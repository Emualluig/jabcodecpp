#include "encoder.h"

#include "reportError.h"

#include <cmath>
#include "ldpc.h"
#include "detector.h"
#include "decoder.h"

#include <optional>
#include "interleave.h"

/**
 * @brief Generate color palettes with more than 8 colors
 * @param color_number the number of colors
 * @param palette the color palette
*/
constexpr static void genColorPalette(jab_int32 color_number, std::unique_ptr<jab_byte[]>& palette) noexcept {
	if (color_number < 8) {
		return;
	}

	jab_int32 vr, vg, vb;	//the number of variable colors for r, g, b channels
	switch(color_number) {
		case 16:
			vr = 4;
			vg = 2;
			vb = 2;
			break;
		case 32:
			vr = 4;
			vg = 4;
			vb = 2;
			break;
		case 64:
			vr = 4;
			vg = 4;
			vb = 4;
			break;
		case 128:
			vr = 8;
			vg = 4;
			vb = 4;
			break;
		case 256:
			vr = 8;
			vg = 8;
			vb = 4;
			break;
		default:
			return;
	}

	//the pixel value interval for r, g, b channels
	jab_float dr = (vr - 1) == 3 ? 85 : 256 / (jab_float)(vr - 1);
	jab_float dg = (vg - 1) == 3 ? 85 : 256 / (jab_float)(vg - 1);
	jab_float db = (vb - 1) == 3 ? 85 : 256 / (jab_float)(vb - 1);

	jab_int32 r, g, b;		//pixel value
	jab_int32 index = 0;	//palette index
	for(jab_int32 i=0; i<vr; i++)
	{
		r = std::min((jab_int32)(dr * i), 255);
		for(jab_int32 j=0; j<vg; j++)
		{
			g = std::min((jab_int32)(dg * j), 255);
			for(jab_int32 k=0; k<vb; k++)
			{
				b = std::min((jab_int32)(db * k), 255);
				palette[index++] = (jab_byte)r;
				palette[index++] = (jab_byte)g;
				palette[index++] = (jab_byte)b;
			}
		}
	}
}

/**
 * @brief Set default color palette
 * @param color_number the number of colors
 * @param palette the color palette
 */
constexpr static void setDefaultPalette(jab_int32 color_number, std::unique_ptr<jab_byte[]>& palette) noexcept {
    if(color_number == 4)
    {
    	memcpy(palette.get() + 0, jab_default_palette.data() + FP0_CORE_COLOR * 3, 3);	//black   000 for 00
    	memcpy(palette.get() + 3, jab_default_palette.data() + 5 * 3, 3);				//magenta 101 for 01
    	memcpy(palette.get() + 6, jab_default_palette.data() + FP2_CORE_COLOR * 3, 3);	//yellow  110 for 10
    	memcpy(palette.get() + 9, jab_default_palette.data() + FP3_CORE_COLOR * 3, 3);	//cyan    011 for 11
    }
    else if(color_number == 8)
    {
        for(jab_int32 i=0; i<color_number*3; i++)
        {
            palette[i] = jab_default_palette[i];
        }
    }
    else
    {
    	genColorPalette(color_number, palette);
    }
}

/**
 * @brief Convert decimal to binary
 * @param dec the decimal value
 * @param bin the data in binary representation
 * @param start_position the position to write in encoded data array
 * @param length the length of the converted binary sequence
 */
static void convert_dec_to_bin(jab_int32 dec, jab_char* bin, jab_int32 start_position, jab_int32 length) noexcept {
	if (dec < 0) {
		dec += 256;
	}
    for (jab_int32 j=0; j < length; j++) {
        jab_char t = dec % 2;
        bin[start_position+length-1-j] = t;
        dec /= 2;
    }
}

/**
 * @brief Create encode object
 * @param color_number the number of module colors
 * @param symbol_number the number of symbols
 * @return the created encode parameter object | NULL: fatal error (out of memory)
 */
Encode createEncode(ColorPaletteType color_number, jab_int32 symbol_number) noexcept {
	if (color_number != ColorPaletteType::FULL && color_number != ColorPaletteType::CMYK) {
		color_number = ColorPaletteType::FULL;
	}

	if (symbol_number < 1 || symbol_number > MAX_SYMBOL_NUMBER) {
		symbol_number = DEFAULT_SYMBOL_NUMBER;
	}

	Encode enc = {
		.color_number = color_number,
		.symbol_number = symbol_number,
		.module_size = DEFAULT_MODULE_SIZE,
		.master_symbol_width = 0,
		.master_symbol_height = 0,
		.palette = std::unique_ptr<jab_byte[]>(new jab_byte[color_number * 3]{}),
		.symbol_versions = std::unique_ptr<jab_vector2d[]>(new jab_vector2d[symbol_number]{}), // allocate memory for symbol versions
		.symbol_ecc_levels = std::unique_ptr<jab_byte[]>(new jab_byte[symbol_number]{}), //set default error correction levels
		.symbol_positions = std::unique_ptr<jab_int32[]>(new jab_int32[symbol_number]{}), // allocate memory for symbol positions
		.symbols = std::unique_ptr<Symbol[]>(new Symbol[symbol_number]), // allocate memory for symbols
	};

    //set default color palette
    setDefaultPalette(enc.color_number, enc.palette);

    return enc;
}

/**
 * @brief Analyze the input data and determine the optimal encoding modes for each character
 * @param input the input character data
 * @param encoded_length the shortest encoding length
 * @return the optimal encoding sequence | NULL: fatal error (out of memory)
 */
static std::optional<std::vector<jab_int32>> analyzeInputData(const std::vector<jab_char>& input, jab_int32& encoded_length) noexcept {
    jab_int32 encode_seq_length = ENC_MAX;
	std::unique_ptr<jab_char[]> seq = std::unique_ptr<jab_char[]>(new jab_char[input.size()]);
	std::unique_ptr<jab_int32[]> curr_seq_len = std::unique_ptr<jab_int32[]>(new jab_int32[(input.size() + 2) * 14]);

	std::unique_ptr<jab_int32[]> prev_mode = std::unique_ptr<jab_int32[]>(new jab_int32[(2 * input.size() + 2) * 14]);
	for (jab_int32 i = 0; i < (2 * input.size() + 2) * 14; i++) {
		prev_mode[i] = ENC_MAX / 2;
	}

	std::array<jab_int32, 28> switch_mode = {};
	std::array<jab_int32, 28> temp_switch_mode = {};
	for (jab_int32 i = 0; i < 28; i++) {
		switch_mode[i] = ENC_MAX / 2;
		temp_switch_mode[i] = ENC_MAX / 2;
	}

    //calculate the shortest encoding sequence
    //initialize start in upper case mode; no previous mode available
    for (jab_int32 k=0;k<7;k++)
    {
        curr_seq_len[k] = curr_seq_len[k+7] = ENC_MAX;
        prev_mode[k] = prev_mode[k+7] = ENC_MAX;
    }

    curr_seq_len[0]=0;
    jab_byte jp_to_nxt_char=0, confirm=0;
    jab_int32 curr_seq_counter=0;
    jab_boolean is_shift=0;
    jab_int32 nb_char=0;
    jab_int32 end_of_loop=input.size();
    jab_int32 prev_mode_index=0;
    for (jab_int32 i=0;i<end_of_loop;i++)
    {
        jab_int32 tmp=input[nb_char];
        jab_int32 tmp1=0;
        if(nb_char+1 < input.size())
            tmp1=input[nb_char+1];
        if(tmp<0)
            tmp=256+tmp;
        if(tmp1<0)
            tmp1=256+tmp1;
        curr_seq_counter++;
        for (jab_int32 j=0;j<JAB_ENCODING_MODES;j++)
        {
            if (jab_enconing_table[tmp][j]>-1 && jab_enconing_table[tmp][j]<64) //check if character is in encoding table
                curr_seq_len[(i+1)*14+j]=curr_seq_len[(i+1)*14+j+7]=character_size[j];
            else if((jab_enconing_table[tmp][j]==-18 && tmp1==10) || (jab_enconing_table[tmp][j]<-18 && tmp1==32))//read next character to decide if encodalbe in current mode
            {
                curr_seq_len[(i+1)*14+j]=curr_seq_len[(i+1)*14+j+7]=character_size[j];
                jp_to_nxt_char=1; //jump to next character
            }
            else //not encodable in this mode
                curr_seq_len[(i+1)*14+j]=curr_seq_len[(i+1)*14+j+7]=ENC_MAX;
        }
        curr_seq_len[(i+1)*14+6]=curr_seq_len[(i+1)*14+13]=character_size[6]; //input sequence can always be encoded by byte mode
        is_shift=0;
        for (jab_int32 j=0;j<14;j++)
        {
            jab_int32 prev=-1;
            jab_int32 len=curr_seq_len[(i+1)*14+j]+curr_seq_len[i*14+j]+latch_shift_to[j][j];
            prev_mode[curr_seq_counter*14+j]=j;
            for (jab_int32 k=0;k<14;k++)
            {
                if((len>=curr_seq_len[(i+1)*14+j]+curr_seq_len[i*14+k]+latch_shift_to[k][j] && k<13) || (k==13 && prev==j))
                {
                    len=curr_seq_len[(i+1)*14+j]+curr_seq_len[i*14+k]+latch_shift_to[k][j];
                    if (temp_switch_mode[2*k]==k)
                        prev_mode[curr_seq_counter*14+j]=temp_switch_mode[2*k+1];
                    else
                        prev_mode[curr_seq_counter*14+j]=k;
                    if (k==13 && prev==j)
                        prev=-1;
                }
            }
            curr_seq_len[(i+1)*14+j]=len;
            //shift back to mode if shift is used
            if (j>6)
            {
                if ((curr_seq_len[(i+1)*14+prev_mode[curr_seq_counter*14+j]]>len ||
                    (jp_to_nxt_char==1 && curr_seq_len[(i+1)*14+prev_mode[curr_seq_counter*14+j]]+character_size[(prev_mode[curr_seq_counter*14+j])%7]>len)) &&
                     j != 13)
                {
                    jab_int32 index=prev_mode[curr_seq_counter*14+j];
                    jab_int32 loop=1;
                    while (index>6 && curr_seq_counter-loop >= 0)
                    {
                        index=prev_mode[(curr_seq_counter-loop)*14+index];
                        loop++;
                    }

                    curr_seq_len[(i+1)*14+index]=len;
                    prev_mode[(curr_seq_counter+1)*14+index]=j;
                    switch_mode[2*index]=index;
                    switch_mode[2*index+1]=j;
                    is_shift=1;
                    if(jp_to_nxt_char==1 && j==11)
                    {
                        confirm=1;
                        prev_mode_index=index;
                    }
                }
                else if ((curr_seq_len[(i+1)*14+prev_mode[curr_seq_counter*14+j]]>len ||
                        (jp_to_nxt_char==1 && curr_seq_len[(i+1)*14+prev_mode[curr_seq_counter*14+j]]+character_size[prev_mode[curr_seq_counter*14+j]%7]>len)) && j == 13 )
                   {
                       curr_seq_len[(i+1)*14+prev_mode[curr_seq_counter*14+j]]=len;
                       prev_mode[(curr_seq_counter+1)*14+prev_mode[curr_seq_counter*14+j]]=j;
                       switch_mode[2*prev_mode[curr_seq_counter*14+j]]=prev_mode[curr_seq_counter*14+j];
                       switch_mode[2*prev_mode[curr_seq_counter*14+j]+1]=j;
                       is_shift=1;
                   }
                if(j!=13)
                    curr_seq_len[(i+1)*14+j]=ENC_MAX;
                else
                    prev=prev_mode[curr_seq_counter*14+j];
            }
        }
        for (jab_int32 j=0;j<28;j++)
        {
            temp_switch_mode[j]=switch_mode[j];
            switch_mode[j]=ENC_MAX/2;
        }

        if(jp_to_nxt_char==1 && confirm==1)
        {
            for (jab_int32 j=0;j<=2*JAB_ENCODING_MODES+1;j++)
            {
                if(j != prev_mode_index)
                    curr_seq_len[(i+1)*14+j]=ENC_MAX;
            }
            nb_char++;
            end_of_loop--;

        }
        jp_to_nxt_char=0;
        confirm=0;
        nb_char++;
    }

    //pick smallest number in last step
    jab_int32 current_mode=0;
    for (jab_int32 j=0;j<=2*JAB_ENCODING_MODES+1;j++)
    {
        if (encode_seq_length>curr_seq_len[(nb_char-(input.size() - end_of_loop)) * 14 + j])
        {
            encode_seq_length=curr_seq_len[(nb_char-(input.size() - end_of_loop)) * 14 + j];
            current_mode=j;
        }
    }
    if(current_mode>6)
        is_shift=1;
    if (is_shift && temp_switch_mode[2*current_mode+1]<14)
        current_mode=temp_switch_mode[2*current_mode+1];

	std::vector<jab_int32> encode_seq = std::vector<jab_int32>(curr_seq_counter + 1 + is_shift);

    //check if byte mode is used more than 15 times in sequence
    //->>length will be increased by 13
    jab_int32 counter=0;
    jab_int32 seq_len=0;
	jab_int32 modeswitch=0;
    encode_seq[curr_seq_counter]=current_mode;//prev_mode[(curr_seq_counter)*14+current_mode];//prev_mode[(curr_seq_counter+is_shift-1)*14+current_mode];
    seq_len+=character_size[encode_seq[curr_seq_counter]%7];
    for (jab_int32 i=curr_seq_counter;i>0;i--)
    {
        if (encode_seq[i]==13 || encode_seq[i]==6)
            counter++;
        else
        {
            if(counter>15)
            {
                encode_seq_length+=13;
                seq_len+=13;

				//--------------------------------
				if(counter>8207) //2^13+15
				{
					if (encode_seq[i]==0 || encode_seq[i]==1 || encode_seq[i]==7 || encode_seq[i]==8)
						modeswitch=11;
					if (encode_seq[i]==2 || encode_seq[i]==9)
						modeswitch=10;
					if (encode_seq[i]==5 || encode_seq[i]==12)
						modeswitch=12;
					jab_int32 remain_in_byte_mode=counter/8207;
					jab_int32 remain_in_byte_mode_residual=counter%8207;
					encode_seq_length+=(remain_in_byte_mode) * modeswitch;
					seq_len+=(remain_in_byte_mode) * modeswitch;
					if(remain_in_byte_mode_residual<16)
					{
						encode_seq_length+=(remain_in_byte_mode-1) * 13;
						seq_len+=(remain_in_byte_mode-1) * 13;
					}
					else
					{
						encode_seq_length+=remain_in_byte_mode * 13;
						seq_len+=remain_in_byte_mode * 13;
					}
					if(remain_in_byte_mode_residual==0)
					{
						encode_seq_length-= modeswitch;
						seq_len-= modeswitch;
					}
				}
				//--------------------------------
				counter=0;
            }
        }
        if (encode_seq[i]<14 && i-1!=0)
        {
            encode_seq[i-1]=prev_mode[i*14+encode_seq[i]];
            seq_len+=character_size[encode_seq[i-1]%7];
            if(encode_seq[i-1]!=encode_seq[i])
                seq_len+=latch_shift_to[encode_seq[i-1]][encode_seq[i]];
        }
        else if (i-1==0)
        {
            encode_seq[i-1]=0;
            if(encode_seq[i-1]!=encode_seq[i])
                seq_len+=latch_shift_to[encode_seq[i-1]][encode_seq[i]];
            if(counter>15)
            {
                encode_seq_length+=13;
                seq_len+=13;

				//--------------------------------
				if(counter>8207) //2^13+15
				{
					modeswitch=11;
					jab_int32 remain_in_byte_mode=counter/8207;
					jab_int32 remain_in_byte_mode_residual=counter%8207;
					encode_seq_length+=remain_in_byte_mode * modeswitch;
					seq_len+=remain_in_byte_mode * modeswitch;
					if(remain_in_byte_mode_residual<16)
					{
						encode_seq_length+=(remain_in_byte_mode-1) * 13;
						seq_len+=(remain_in_byte_mode-1) * 13;
					}
					else
					{
						encode_seq_length+=remain_in_byte_mode * 13;
						seq_len+=remain_in_byte_mode * 13;
					}
					if(remain_in_byte_mode_residual==0)
					{
						encode_seq_length-= modeswitch;
						seq_len-= modeswitch;
					}
				}
				//--------------------------------
				counter=0;
            }
        }
		else {
			return std::nullopt;
		}
    }
    encoded_length = encode_seq_length;
    return encode_seq;
}

/**
 * @brief Check if master symbol shall be encoded in default mode
 * @param enc the encode parameters
 * @return JAB_SUCCESS | JAB_FAILURE
*/
static jab_boolean isDefaultMode(const Encode& enc) noexcept {
	if(enc.color_number == 8 && (enc.symbol_ecc_levels[0] == 0 || enc.symbol_ecc_levels[0] == DEFAULT_ECC_LEVEL))
	{
		return JAB_SUCCESS;
	}
	return JAB_FAILURE;
}

/**
 * @brief Calculate the (encoded) metadata length
 * @param enc the encode parameters
 * @param index the symbol index
 * @return the metadata length (encoded length for master symbol)
*/
jab_int32 getMetadataLength(const Encode& enc, jab_int32 index)
{
    jab_int32 length = 0;

    if (index == 0) //master symbol, the encoded length
    {
    	//default mode, no metadata
    	if(isDefaultMode(enc))
		{
			length = 0;
		}
		else
		{
			//Part I
			length += MASTER_METADATA_PART1_LENGTH;
			//Part II
			length += MASTER_METADATA_PART2_LENGTH;
		}
    }
    else //slave symbol, the original net length
    {
    	//Part I
        length += 2;
        //Part II
        jab_int32 host_index = enc.symbols[index].host;
        //V in Part II, compare symbol shape and size with host symbol
        if (enc.symbol_versions[index].x != enc.symbol_versions[host_index].x || enc.symbol_versions[index].y != enc.symbol_versions[host_index].y)
		{
			length += 5;
		}
        //E in Part II
        if (enc.symbol_ecc_levels[index] != enc.symbol_ecc_levels[host_index])
        {
            length += 6;
        }
    }
    return length;
}

/**
 * @brief Calculate the data capacity of a symbol
 * @param enc the encode parameters
 * @param index the symbol index
 * @return the data capacity
 */
static jab_int32 getSymbolCapacity(const Encode& enc, jab_int32 index) noexcept {
	//number of modules for finder patterns
    jab_int32 nb_modules_fp;
    if(index == 0)	//master symbol
	{
		nb_modules_fp = 4 * 17;
	}
	else			//slave symbol
	{
		nb_modules_fp = 4 * 7;
	}
    //number of modules for color palette
    jab_int32 nb_modules_palette = enc.color_number > 64 ? (64-2)*COLOR_PALETTE_NUMBER : (enc.color_number-2)*COLOR_PALETTE_NUMBER;
	//number of modules for alignment pattern
	jab_int32 side_size_x = VERSION2SIZE(enc.symbol_versions[index].x);
	jab_int32 side_size_y = VERSION2SIZE(enc.symbol_versions[index].y);
	jab_int32 number_of_aps_x = jab_ap_num[enc.symbol_versions[index].x - 1];
	jab_int32 number_of_aps_y = jab_ap_num[enc.symbol_versions[index].y - 1];
	jab_int32 nb_modules_ap = (number_of_aps_x * number_of_aps_y - 4) * 7;
	//number of modules for metadata
	jab_int32 nb_of_bpm = std::log2(enc.color_number);
	jab_int32 nb_modules_metadata = 0;
	if(index == 0)	//master symbol
	{
		jab_int32 nb_metadata_bits = getMetadataLength(enc, index);
		if(nb_metadata_bits > 0)
		{
			nb_modules_metadata = (nb_metadata_bits - MASTER_METADATA_PART1_LENGTH) / nb_of_bpm; //only modules for PartII
			if((nb_metadata_bits - MASTER_METADATA_PART1_LENGTH) % nb_of_bpm != 0)
			{
				nb_modules_metadata++;
			}
			nb_modules_metadata += MASTER_METADATA_PART1_MODULE_NUMBER; //add modules for PartI
		}
	}
	jab_int32 capacity = (side_size_x*side_size_y - nb_modules_fp - nb_modules_ap - nb_modules_palette - nb_modules_metadata) * nb_of_bpm;
	return capacity;
}

/**
 * @brief Get the optimal error correction capability
 * @param capacity the symbol capacity
 * @param net_data_length the original data length
 * @param wcwr the LPDC parameters wc and wr
 * @return JAB_SUCCESS | JAB_FAILURE
 */
void getOptimalECC(jab_int32 capacity, jab_int32 net_data_length, jab_int32* wcwr)
{
	jab_float min = capacity;
	for (jab_int32 k=3; k<=6+2; k++)
	{
		for (jab_int32 j=k+1; j<=6+3; j++)
		{
			jab_int32 dist = (capacity/j)*j - (capacity/j)*k - net_data_length; //max_gross_payload = floor(capacity / wr) * wr
			if(dist<min && dist>=0)
			{
				wcwr[1] = j;
				wcwr[0] = k;
				min = dist;
			}
		}
	}
}

/**
 * @brief Encode the input data
 * @param data the character input data
 * @param encoded_length the optimal encoding length
 * @param encode_seq the optimal encoding sequence
 * @return the encoded data | NULL if failed
 */
static std::optional<std::vector<jab_char>> encodeData(const std::vector<jab_char>& data, jab_int32 encoded_length, std::vector<jab_int32>& encode_seq) noexcept {
	std::vector<jab_char> encoded_data = std::vector<jab_char>(encoded_length);

    jab_int32 counter=0;
    jab_boolean shift_back=0;
    jab_int32 position=0;
    jab_int32 current_encoded_length=0;
    jab_int32 end_of_loop=data.size();
    jab_int32 byte_offset=0;
    jab_int32 byte_counter=0;
	jab_int32 factor=1;
    //encoding starts in upper case mode
    for (jab_int32 i=0;i<end_of_loop;i++)
    {
        jab_int32 tmp=data[current_encoded_length];
        if (tmp < 0)
            tmp+=256;
        if (position<encoded_length)
        {
            jab_int32 decimal_value;
            //check if mode is switched
            if (encode_seq[counter] != encode_seq[counter+1])
            {
                //encode mode switch
                jab_int32 length=latch_shift_to[encode_seq[counter]][encode_seq[counter+1]];
                if(encode_seq[counter+1] == 6 || encode_seq[counter+1] == 13)
                    length-=4;
                if(length < ENC_MAX)
                    convert_dec_to_bin(mode_switch[encode_seq[counter]][encode_seq[counter+1]],encoded_data.data(), position, length);
                else
                {
                    reportError("Encoding data failed");
                    return std::nullopt;
                }
                position+=latch_shift_to[encode_seq[counter]][encode_seq[counter+1]];
                if(encode_seq[counter+1] == 6 || encode_seq[counter+1] == 13)
                    position-=4;
                //check if latch or shift
                if((encode_seq[counter+1]>6 && encode_seq[counter+1]<=13) || (encode_seq[counter+1]==13 && encode_seq[counter+2]!=13))
                    shift_back=1;//remember to shift back to mode from which was invoked
            }
            //if not byte mode
            if (encode_seq[counter+1]%7 != 6)// || end_of_loop-1 == i)
            {
                if(jab_enconing_table[tmp][encode_seq[counter+1]%7]>-1 && character_size[encode_seq[counter+1]%7] < ENC_MAX)
                {
                    //encode character
                    convert_dec_to_bin(jab_enconing_table[tmp][encode_seq[counter+1]%7],encoded_data.data(), position, character_size[encode_seq[counter + 1] % 7]);
                    position+=character_size[encode_seq[counter+1]%7];
                    counter++;
                }
                else if (jab_enconing_table[tmp][encode_seq[counter+1]%7]<-1)
                {
                    jab_int32 tmp1=data[current_encoded_length+1];
                    if (tmp1 < 0)
                        tmp1+=256;
                    //read next character to see if more efficient encoding possible
                    if (((tmp==44 || tmp== 46 || tmp==58) && tmp1==32) || (tmp==13 && tmp1==10))
                        decimal_value=abs(jab_enconing_table[tmp][encode_seq[counter+1]%7]);
                    else if (tmp==13 && tmp1!=10)
                        decimal_value=18;
                    else
                    {
                        reportError("Encoding data failed");
                        return std::nullopt;
                    }
                    if (character_size[encode_seq[counter+1]%7] < ENC_MAX)
                    convert_dec_to_bin(decimal_value,encoded_data.data(), position, character_size[encode_seq[counter + 1] % 7]);
                    position+=character_size[encode_seq[counter+1]%7];
                    counter++;
                    end_of_loop--;
                    current_encoded_length++;
                }
                else
                {
                    reportError("Encoding data failed");
                    return std::nullopt;
                }
            }
            else
            {
                //byte mode
                if(encode_seq[counter] != encode_seq[counter+1])
                {
                    //loop over sequence to check how many characters in byte mode follow
                    byte_counter=0;
                    for(jab_int32 byte_loop=counter+1;byte_loop<=end_of_loop;byte_loop++)
                    {
                        if(encode_seq[byte_loop]==6 || encode_seq[byte_loop]==13)
                            byte_counter++;
                        else
                            break;
                    }
                    convert_dec_to_bin(byte_counter > 15 ? 0 : byte_counter,encoded_data.data(), position, 4);
                    position+=4;
                    if(byte_counter > 15)
                    {
						if(byte_counter <= 8207)//8207=2^13+15; if number of bytes exceeds 8207, encoder shall shift to byte mode again from upper case mode && byte_counter < 8207
						{
							convert_dec_to_bin(byte_counter-15-1,encoded_data.data(), position, 13);
						}
						else
						{
							convert_dec_to_bin(8191,encoded_data.data(), position, 13);
						}
                        position+=13;
                    }
                    byte_offset=byte_counter;
                }
				if(byte_offset-byte_counter==factor*8207) //byte mode exceeds 2^13 + 15
				{
					if(encode_seq[counter-(byte_offset-byte_counter)]==0 || encode_seq[counter-(byte_offset-byte_counter)]==7 || encode_seq[counter-(byte_offset-byte_counter)]==1|| encode_seq[counter-(byte_offset-byte_counter)]==8)
					{
						convert_dec_to_bin(124,encoded_data.data(), position, 7);// shift from upper case to byte
						position+=7;
					}
					if(encode_seq[counter-(byte_offset-byte_counter)]==2 || encode_seq[counter-(byte_offset-byte_counter)]==9)
					{
						convert_dec_to_bin(60,encoded_data.data(), position, 5);// shift from numeric to byte
						position+=5;
					}
					if(encode_seq[counter-(byte_offset-byte_counter)]==5 || encode_seq[counter-(byte_offset-byte_counter)]==12)
					{
						convert_dec_to_bin(252,encoded_data.data(), position, 8);// shift from alphanumeric to byte
						position+=8;
					}
					convert_dec_to_bin(byte_counter > 15 ? 0 : byte_counter,encoded_data.data(), position, 4); //write the first 4 bits
					position+=4;
					if(byte_counter > 15) //if more than 15 bytes -> use the next 13 bits to wirte the length
					{
						if(byte_counter <= 8207)//8207=2^13+15; if number of bytes exceeds 8207, encoder shall shift to byte mode again from upper case mode && byte_counter < 8207
						{
							convert_dec_to_bin(byte_counter-15-1,encoded_data.data(), position, 13);
						}
						else //number exceeds 2^13 + 15
						{
							convert_dec_to_bin(8191,encoded_data.data(), position, 13);
						}
						position+=13;
					}
					factor++;
				}
                if (character_size[encode_seq[counter+1]%7] < ENC_MAX)
                    convert_dec_to_bin(tmp,encoded_data.data(), position, character_size[encode_seq[counter + 1] % 7]);
                else
                {
                    reportError("Encoding data failed");
                    return std::nullopt;
                }
                position+=character_size[encode_seq[counter+1]%7];
                counter++;
                byte_counter--;
            }
            //shift back to mode from which mode was invoked
            if (shift_back && byte_counter==0)
            {
                if(byte_offset==0)
                    encode_seq[counter]=encode_seq[counter-1];
                else
                    encode_seq[counter]=encode_seq[counter-byte_offset];
                shift_back=0;
                byte_offset=0;
            }

        }
        else
        {
            reportError("Encoding data failed");
            return std::nullopt;
        }
        current_encoded_length++;
    }
    return encoded_data;
}

/**
 * @brief Encode metadata
 * @param enc the encode parameters
 * @return JAB_SUCCESS | JAB_FAILURE
*/
static jab_boolean encodeMasterMetadata(Encode& enc) noexcept {
	jab_int32 partI_length 	= MASTER_METADATA_PART1_LENGTH/2;	//partI net length
	jab_int32 partII_length	= MASTER_METADATA_PART2_LENGTH/2;	//partII net length
	jab_int32 V_length = 10;
	jab_int32 E_length = 6;
	jab_int32 MSK_length = 3;
	//set master metadata variables
	jab_int32 Nc = std::log2(enc.color_number) - 1;
	jab_int32 V = ((enc.symbol_versions[0].x -1) << 5) + (enc.symbol_versions[0].y - 1);
	jab_int32 E1 = enc.symbols[0].wcwr[0] - 3;
	jab_int32 E2 = enc.symbols[0].wcwr[1] - 4;
	jab_int32 MSK = DEFAULT_MASKING_REFERENCE;

	//write each part of master metadata
	//Part I
	jab_char* partI = new jab_char[partI_length];
	convert_dec_to_bin(Nc, partI, 0, partI_length);

	//Part II
	jab_char* partII = new jab_char[partII_length];
	convert_dec_to_bin(V,   partII, 0, V_length);
	convert_dec_to_bin(E1,  partII, V_length, 3);
	convert_dec_to_bin(E2,  partII, V_length+3, 3);
	convert_dec_to_bin(MSK, partII, V_length+E_length, MSK_length);

	//encode each part of master metadata
	jab_int32 wcwr[2] = {2, -1};

	//Part I
	std::vector<jab_char> encoded_partI   = encodeLDPC(partI, wcwr, partI_length); // FIX1: Doubles the length since wr < 0
	delete[] partI;

	//Part II
	std::vector<jab_char> encoded_partII  = encodeLDPC(partII, wcwr, partII_length); // FIX1: Doubles the length since wr < 0
	delete[] partII;

	jab_int32 encoded_metadata_length = encoded_partI.size() + encoded_partII.size(); // FIX1: 2 * 3 + 2 * 19 = 44
	enc.symbols[0].metadata = std::vector<jab_char>(encoded_metadata_length);

	//copy encoded parts into metadata
	memcpy(enc.symbols[0].metadata.data(), encoded_partI.data(), encoded_partI.size());
	memcpy(enc.symbols[0].metadata.data() + encoded_partI.size(), encoded_partII.data(), encoded_partII.size());

	// FIX1: The master symbol's length has been set to 44
    return JAB_SUCCESS;
}

#include <array>
#include "mask.h"
/**
 * @brief Update master symbol metadata PartII if the default masking reference is changed
 * @param enc the encode parameter
 * @param mask_ref the masking reference
 * @return JAB_SUCCESS | JAB_FAILURE
*/
static jab_boolean updateMasterMetadataPartII(Encode& enc, jab_int32 mask_ref) noexcept {
	std::array<jab_char, MASTER_METADATA_PART2_LENGTH / 2> partII = {};

	//set V and E
	jab_int32 V_length = 10;
	jab_int32 E_length = 6;
	jab_int32 MSK_length = 3;
	jab_int32 V = ((enc.symbol_versions[0].x -1) << 5) + (enc.symbol_versions[0].y - 1);
	jab_int32 E1 = enc.symbols[0].wcwr[0] - 3;
	jab_int32 E2 = enc.symbols[0].wcwr[1] - 4;
	convert_dec_to_bin(V,   partII.data(), 0, V_length);
	convert_dec_to_bin(E1,  partII.data(), V_length, 3);
	convert_dec_to_bin(E2,  partII.data(), V_length+3, 3);

	//update masking reference in PartII
	convert_dec_to_bin(mask_ref, partII.data(), V_length+E_length, MSK_length);

	//encode new PartII
	jab_int32 wcwr[2] = {2, -1};
	std::vector<jab_char> encoded_partII = encodeLDPC(partII.data(), wcwr, partII.size());

	//update metadata
	memcpy(enc.symbols[0].metadata.data() + MASTER_METADATA_PART1_LENGTH, encoded_partII.data(), encoded_partII.size());

	return JAB_SUCCESS;
}

/**
 * @brief Update master symbol metadata PartII if the default masking reference is changed
 * @param enc the encode parameter
*/
static void placeMasterMetadataPartII(Encode& enc) noexcept {
    //rewrite metadata in master with mask information
	jab_int32 nb_of_bits_per_mod = std::log2(enc.color_number);
    jab_int32 x = MASTER_METADATA_X;
    jab_int32 y = MASTER_METADATA_Y;
    jab_int32 module_count = 0;
    //skip PartI and color palette
    jab_int32 color_palette_size = std::min(enc.color_number-2, 64-2);
    jab_int32 module_offset = MASTER_METADATA_PART1_MODULE_NUMBER + color_palette_size*COLOR_PALETTE_NUMBER;
    for(jab_int32 i=0; i<module_offset; i++)
	{
		module_count++;
        getNextMetadataModuleInMaster(enc.symbols[0].side_size.y, enc.symbols[0].side_size.x, module_count, x, y);
	}
	//update PartII
	jab_int32 partII_bit_start = MASTER_METADATA_PART1_LENGTH;
	jab_int32 partII_bit_end = MASTER_METADATA_PART1_LENGTH + MASTER_METADATA_PART2_LENGTH;
	jab_int32 metadata_index = partII_bit_start;

	// FIX1: This function is called after encodeMasterMetadata, which set the master symbol's length to 44.
	//       But partII_bit_end is 44, which means there was a array overrun. This was fixed by replacing <= with <
	while(metadata_index < partII_bit_end) { // FIX1: changed <= to <
    	jab_byte color_index = enc.symbols[0].matrix[y*enc.symbols[0].side_size.x + x];
		for(jab_int32 j = 0; j< nb_of_bits_per_mod; j++)
		{
			if(metadata_index < partII_bit_end) // FIX1: also changed <= to < here
			{
				// FIX1: Here is the potential buffer over run. 
				jab_byte bit = enc.symbols[0].metadata[metadata_index];

				if(bit == 0)
					color_index &= ~(1UL << (nb_of_bits_per_mod-1-j));
				else
					color_index |= 1UL << (nb_of_bits_per_mod-1-j);
				metadata_index++;
			}
			else
				break;
		}
        enc.symbols[0].matrix[y*enc.symbols[0].side_size.x + x] = color_index;
        module_count++;
        getNextMetadataModuleInMaster(enc.symbols[0].side_size.y, enc.symbols[0].side_size.x, module_count, x, y);
    }
}

/**
 * @brief Get color index for the color palette
 * @param index the color index in the palette
 * @param index_size the size of index
 * @param color_number the number of colors
*/
static void getColorPaletteIndex(std::vector<jab_byte>& index, jab_int32 index_size, jab_int32 color_number) noexcept {
	for(jab_int32 i = 0; i< index_size; i++) {
		index[i] = i;
	}

	if (color_number != 128 && color_number != 256) {
		return;
	}

	std::unique_ptr<jab_byte[]> tmpVector = std::unique_ptr<jab_byte[]>(new jab_byte[color_number]);
	for(jab_int32 i=0; i< color_number; i++) {
		tmpVector[i] = i;
	}
	jab_byte* tmp = tmpVector.get();

	jab_byte* ptr = index.data();
	if(color_number == 128) {
		memcpy(ptr + 0,  tmp + 0, 16);
		memcpy(ptr + 16, tmp + 32, 16);
		memcpy(ptr + 32, tmp + 80, 16);
		memcpy(ptr + 48, tmp + 112, 16);
	}
	else if(color_number == 256) {
		memcpy(ptr + 0, tmp + 0,  4);
		memcpy(ptr + 4, tmp + 8,  4);
		memcpy(ptr + 8, tmp + 20, 4);
		memcpy(ptr + 12,tmp + 28, 4);

		memcpy(ptr + 16, tmp + 64, 4);
		memcpy(ptr + 20, tmp + 72, 4);
		memcpy(ptr + 24, tmp + 84, 4);
		memcpy(ptr + 28, tmp + 92, 4);

		memcpy(ptr + 32, tmp + 160, 4);
		memcpy(ptr + 36, tmp + 168, 4);
		memcpy(ptr + 40, tmp + 180, 4);
		memcpy(ptr + 44, tmp + 188, 4);

		memcpy(ptr + 48, tmp + 224, 4);
		memcpy(ptr + 52, tmp + 232, 4);
		memcpy(ptr + 56, tmp + 244, 4);
		memcpy(ptr + 60, tmp + 252, 4);
	}
}

/**
 * @brief Create symbol matrix
 * @param enc the encode parameter
 * @param index the symbol index
 * @param ecc_encoded_data encoded data
 * @return JAB_SUCCESS | JAB_FAILURE
*/
static jab_boolean createMatrix(Encode& enc, jab_int32 index, std::vector<jab_char>& ecc_encoded_data) noexcept {
    //Allocate matrix
	enc.symbols[index].matrix = std::unique_ptr<jab_byte[]>(new jab_byte[enc.symbols[index].side_size.x * enc.symbols[index].side_size.y]{});

    //Allocate boolean matrix
    enc.symbols[index].data_map = std::vector<jab_byte>(enc.symbols[index].side_size.x * enc.symbols[index].side_size.y, 1);

    //set alignment patterns
    jab_int32 Nc = std::log2(enc.color_number) - 1;
	jab_byte apx_core_color = apx_core_color_index[Nc];
	jab_byte apx_peri_color = apn_core_color_index[Nc];
	jab_int32 side_ver_x_index = SIZE2VERSION(enc.symbols[index].side_size.x) - 1;
	jab_int32 side_ver_y_index = SIZE2VERSION(enc.symbols[index].side_size.y) - 1;
    for(jab_int32 x=0; x<jab_ap_num[side_ver_x_index]; x++)
    {
    	jab_byte left;
        if (x%2 == 1)
            left=0;
        else
            left=1;
        for(jab_int32 y=0; y<jab_ap_num[side_ver_y_index]; y++)
        {
            jab_int32 x_offset = jab_ap_pos[side_ver_x_index][x] - 1;
            jab_int32 y_offset = jab_ap_pos[side_ver_y_index][y] - 1;
            //left alignment patterns
            if(	left == 1
				&& (x != 0 || y != 0)
				&& (x != 0 || y != jab_ap_num[side_ver_y_index]-1)
				&& (x != jab_ap_num[side_ver_x_index]-1 || y != 0)
				&& (x != jab_ap_num[side_ver_x_index]-1 || y != jab_ap_num[side_ver_y_index]-1))
            {
            	enc.symbols[index].matrix[(y_offset-1)*enc.symbols[index].side_size.x + x_offset-1]=
				enc.symbols[index].matrix[(y_offset-1)*enc.symbols[index].side_size.x + x_offset  ]=
				enc.symbols[index].matrix[(y_offset  )*enc.symbols[index].side_size.x + x_offset-1]=
				enc.symbols[index].matrix[(y_offset  )*enc.symbols[index].side_size.x + x_offset+1]=
				enc.symbols[index].matrix[(y_offset+1)*enc.symbols[index].side_size.x + x_offset  ]=
				enc.symbols[index].matrix[(y_offset+1)*enc.symbols[index].side_size.x + x_offset+1]=apx_peri_color;
				enc.symbols[index].matrix[(y_offset  )*enc.symbols[index].side_size.x + x_offset  ]=apx_core_color;

				enc.symbols[index].data_map[(y_offset-1)*enc.symbols[index].side_size.x + x_offset-1]=
				enc.symbols[index].data_map[(y_offset-1)*enc.symbols[index].side_size.x + x_offset  ]=
				enc.symbols[index].data_map[(y_offset  )*enc.symbols[index].side_size.x + x_offset-1]=
				enc.symbols[index].data_map[(y_offset  )*enc.symbols[index].side_size.x + x_offset+1]=
				enc.symbols[index].data_map[(y_offset+1)*enc.symbols[index].side_size.x + x_offset  ]=
				enc.symbols[index].data_map[(y_offset+1)*enc.symbols[index].side_size.x + x_offset+1]=
				enc.symbols[index].data_map[(y_offset  )*enc.symbols[index].side_size.x + x_offset  ]=0;
            }
            //right alignment patterns
            else if(left == 0
					&& (x != 0 || y != 0)
					&& (x != 0 || y != jab_ap_num[side_ver_y_index]-1)
					&& (x != jab_ap_num[side_ver_x_index]-1 || y != 0)
					&& (x != jab_ap_num[side_ver_x_index]-1 || y != jab_ap_num[side_ver_y_index]-1))
            {
            	enc.symbols[index].matrix[(y_offset-1)*enc.symbols[index].side_size.x + x_offset+1]=
				enc.symbols[index].matrix[(y_offset-1)*enc.symbols[index].side_size.x + x_offset  ]=
				enc.symbols[index].matrix[(y_offset  )*enc.symbols[index].side_size.x + x_offset-1]=
				enc.symbols[index].matrix[(y_offset  )*enc.symbols[index].side_size.x + x_offset+1]=
				enc.symbols[index].matrix[(y_offset+1)*enc.symbols[index].side_size.x + x_offset  ]=
				enc.symbols[index].matrix[(y_offset+1)*enc.symbols[index].side_size.x + x_offset-1]=apx_peri_color;
				enc.symbols[index].matrix[(y_offset  )*enc.symbols[index].side_size.x + x_offset  ]=apx_core_color;

				enc.symbols[index].data_map[(y_offset-1)*enc.symbols[index].side_size.x + x_offset+1]=
				enc.symbols[index].data_map[(y_offset-1)*enc.symbols[index].side_size.x + x_offset  ]=
				enc.symbols[index].data_map[(y_offset  )*enc.symbols[index].side_size.x + x_offset-1]=
				enc.symbols[index].data_map[(y_offset  )*enc.symbols[index].side_size.x + x_offset+1]=
				enc.symbols[index].data_map[(y_offset+1)*enc.symbols[index].side_size.x + x_offset  ]=
				enc.symbols[index].data_map[(y_offset+1)*enc.symbols[index].side_size.x + x_offset-1]=
				enc.symbols[index].data_map[(y_offset  )*enc.symbols[index].side_size.x + x_offset  ]=0;
            }
            if (left==0)
                left=1;
            else
                left=0;
        }
    }

    //outer layers of finder pattern for master symbol
    if(index == 0)
    {
        //if k=0 center, k=1 first layer, k=2 second layer
        for(jab_int32 k=0;k<3;k++)
        {
            for(jab_int32 i=0;i<k+1;i++)
            {
                for(jab_int32 j=0;j<k+1;j++)
                {
                    if (i==k || j==k)
                    {
                        jab_byte fp0_color_index, fp1_color_index, fp2_color_index, fp3_color_index;
						fp0_color_index = (k%2) ? fp3_core_color_index[Nc] : fp0_core_color_index[Nc];
						fp1_color_index = (k%2) ? fp2_core_color_index[Nc] : fp1_core_color_index[Nc];
						fp2_color_index = (k%2) ? fp1_core_color_index[Nc] : fp2_core_color_index[Nc];
						fp3_color_index = (k%2) ? fp0_core_color_index[Nc] : fp3_core_color_index[Nc];

						//upper pattern
                        enc.symbols[index].matrix[(DISTANCE_TO_BORDER-(i+1))*enc.symbols[index].side_size.x+DISTANCE_TO_BORDER-j-1]=
						enc.symbols[index].matrix[(DISTANCE_TO_BORDER+(i-1))*enc.symbols[index].side_size.x+DISTANCE_TO_BORDER+j-1]=fp0_color_index;
                        enc.symbols[index].data_map[(DISTANCE_TO_BORDER-(i+1))*enc.symbols[index].side_size.x+DISTANCE_TO_BORDER-j-1]=
                        enc.symbols[index].data_map[(DISTANCE_TO_BORDER+(i-1))*enc.symbols[index].side_size.x+DISTANCE_TO_BORDER+j-1]=0;

                        enc.symbols[index].matrix[(DISTANCE_TO_BORDER-(i+1))*enc.symbols[index].side_size.x+enc.symbols[index].side_size.x-(DISTANCE_TO_BORDER-1)-j-1]=
                        enc.symbols[index].matrix[(DISTANCE_TO_BORDER+(i-1))*enc.symbols[index].side_size.x+enc.symbols[index].side_size.x-(DISTANCE_TO_BORDER-1)+j-1]=fp1_color_index;
                        enc.symbols[index].data_map[(DISTANCE_TO_BORDER-(i+1))*enc.symbols[index].side_size.x+enc.symbols[index].side_size.x-(DISTANCE_TO_BORDER-1)-j-1]=
                        enc.symbols[index].data_map[(DISTANCE_TO_BORDER+(i-1))*enc.symbols[index].side_size.x+enc.symbols[index].side_size.x-(DISTANCE_TO_BORDER-1)+j-1]=0;

                        //lower pattern
                        enc.symbols[index].matrix[(enc.symbols[index].side_size.y-DISTANCE_TO_BORDER+i)*enc.symbols[index].side_size.x+enc.symbols[index].side_size.x-(DISTANCE_TO_BORDER-1)-j-1]=
                        enc.symbols[index].matrix[(enc.symbols[index].side_size.y-DISTANCE_TO_BORDER-i)*enc.symbols[index].side_size.x+enc.symbols[index].side_size.x-(DISTANCE_TO_BORDER-1)+j-1]=fp2_color_index;
                        enc.symbols[index].data_map[(enc.symbols[index].side_size.y-DISTANCE_TO_BORDER+i)*enc.symbols[index].side_size.x+enc.symbols[index].side_size.x-(DISTANCE_TO_BORDER-1)-j-1]=
                        enc.symbols[index].data_map[(enc.symbols[index].side_size.y-DISTANCE_TO_BORDER-i)*enc.symbols[index].side_size.x+enc.symbols[index].side_size.x-(DISTANCE_TO_BORDER-1)+j-1]=0;

                        enc.symbols[index].matrix[(enc.symbols[index].side_size.y-DISTANCE_TO_BORDER+i)*enc.symbols[index].side_size.x+(DISTANCE_TO_BORDER)-j-1]=
                        enc.symbols[index].matrix[(enc.symbols[index].side_size.y-DISTANCE_TO_BORDER-i)*enc.symbols[index].side_size.x+(DISTANCE_TO_BORDER)+j-1]=fp3_color_index;
                        enc.symbols[index].data_map[(enc.symbols[index].side_size.y-DISTANCE_TO_BORDER+i)*enc.symbols[index].side_size.x+(DISTANCE_TO_BORDER)-j-1]=
                        enc.symbols[index].data_map[(enc.symbols[index].side_size.y-DISTANCE_TO_BORDER-i)*enc.symbols[index].side_size.x+(DISTANCE_TO_BORDER)+j-1]=0;
                    }
                }
            }
        }
    }
    else //finder alignments in slave
    {
        //if k=0 center, k=1 first layer
        for(jab_int32 k=0;k<2;k++)
        {
            for(jab_int32 i=0;i<k+1;i++)
            {
                for(jab_int32 j=0;j<k+1;j++)
                {
                    if (i==k || j==k)
                    {
                    	jab_byte ap0_color_index, ap1_color_index, ap2_color_index, ap3_color_index;
                        ap0_color_index =
						ap1_color_index =
						ap2_color_index =
						ap3_color_index = (k%2) ? apx_core_color_index[Nc] : apn_core_color_index[Nc];
                        //upper pattern
                        enc.symbols[index].matrix[(DISTANCE_TO_BORDER-(i+1))*enc.symbols[index].side_size.x+DISTANCE_TO_BORDER-j-1]=
                        enc.symbols[index].matrix[(DISTANCE_TO_BORDER+(i-1))*enc.symbols[index].side_size.x+DISTANCE_TO_BORDER+j-1]=ap0_color_index;
                        enc.symbols[index].data_map[(DISTANCE_TO_BORDER-(i+1))*enc.symbols[index].side_size.x+DISTANCE_TO_BORDER-j-1]=
                        enc.symbols[index].data_map[(DISTANCE_TO_BORDER+(i-1))*enc.symbols[index].side_size.x+DISTANCE_TO_BORDER+j-1]=0;

                        enc.symbols[index].matrix[(DISTANCE_TO_BORDER-(i+1))*enc.symbols[index].side_size.x+enc.symbols[index].side_size.x-(DISTANCE_TO_BORDER-1)-j-1]=
                        enc.symbols[index].matrix[(DISTANCE_TO_BORDER+(i-1))*enc.symbols[index].side_size.x+enc.symbols[index].side_size.x-(DISTANCE_TO_BORDER-1)+j-1]=ap1_color_index;
                        enc.symbols[index].data_map[(DISTANCE_TO_BORDER-(i+1))*enc.symbols[index].side_size.x+enc.symbols[index].side_size.x-(DISTANCE_TO_BORDER-1)-j-1]=
                        enc.symbols[index].data_map[(DISTANCE_TO_BORDER+(i-1))*enc.symbols[index].side_size.x+enc.symbols[index].side_size.x-(DISTANCE_TO_BORDER-1)+j-1]=0;

                        //lower pattern
                        enc.symbols[index].matrix[(enc.symbols[index].side_size.y-DISTANCE_TO_BORDER+i)*enc.symbols[index].side_size.x+enc.symbols[index].side_size.x-(DISTANCE_TO_BORDER-1)-j-1]=
                        enc.symbols[index].matrix[(enc.symbols[index].side_size.y-DISTANCE_TO_BORDER-i)*enc.symbols[index].side_size.x+enc.symbols[index].side_size.x-(DISTANCE_TO_BORDER-1)+j-1]=ap2_color_index;
                        enc.symbols[index].data_map[(enc.symbols[index].side_size.y-DISTANCE_TO_BORDER+i)*enc.symbols[index].side_size.x+enc.symbols[index].side_size.x-(DISTANCE_TO_BORDER-1)-j-1]=
						enc.symbols[index].data_map[(enc.symbols[index].side_size.y-DISTANCE_TO_BORDER-i)*enc.symbols[index].side_size.x+enc.symbols[index].side_size.x-(DISTANCE_TO_BORDER-1)+j-1]=0;

                        enc.symbols[index].matrix[(enc.symbols[index].side_size.y-DISTANCE_TO_BORDER+i)*enc.symbols[index].side_size.x+(DISTANCE_TO_BORDER)-j-1]=
                        enc.symbols[index].matrix[(enc.symbols[index].side_size.y-DISTANCE_TO_BORDER-i)*enc.symbols[index].side_size.x+(DISTANCE_TO_BORDER)+j-1]=ap3_color_index;
                        enc.symbols[index].data_map[(enc.symbols[index].side_size.y-DISTANCE_TO_BORDER+i)*enc.symbols[index].side_size.x+(DISTANCE_TO_BORDER)-j-1]=
                        enc.symbols[index].data_map[(enc.symbols[index].side_size.y-DISTANCE_TO_BORDER-i)*enc.symbols[index].side_size.x+(DISTANCE_TO_BORDER)+j-1]=0;
                    }
                }
            }
        }
    }

    //Metadata and color palette placement
    jab_int32 nb_of_bits_per_mod = std::log2(enc.color_number);
    jab_int32 color_index;
    jab_int32 module_count = 0;
    jab_int32 x;
    jab_int32 y;

    //get color index for color palette
    jab_int32 palette_index_size = enc.color_number > 64 ? 64 : enc.color_number;
	std::vector<jab_byte> palette_index = std::vector<jab_byte>(palette_index_size);

	getColorPaletteIndex(palette_index, palette_index_size, enc.color_number);

    if(index == 0)//place metadata and color palette in master symbol
    {
        x = MASTER_METADATA_X;
        y = MASTER_METADATA_Y;
		int metadata_index = 0;
        //metadata Part I
        if(!isDefaultMode(enc))
		{
			while(metadata_index < enc.symbols[index].metadata.size() && metadata_index < MASTER_METADATA_PART1_LENGTH)
			{
				//Read 3 bits from encoded PartI each time
				jab_byte bit1 = enc.symbols[index].metadata[metadata_index + 0];
				jab_byte bit2 = enc.symbols[index].metadata[metadata_index + 1];
				jab_byte bit3 = enc.symbols[index].metadata[metadata_index + 2];
				jab_int32 val = (bit1 << 2) + (bit2 << 1) + bit3;
				//place two modules according to the value of every 3 bits
                for(jab_int32 i=0; i<2; i++)
				{
					color_index = nc_color_encode_table[val][i] % enc.color_number;
					enc.symbols[index].matrix  [y*enc.symbols[index].side_size.x + x] = color_index;
					enc.symbols[index].data_map[y*enc.symbols[index].side_size.x + x] = 0;
					module_count++;
					getNextMetadataModuleInMaster(enc.symbols[index].side_size.y, enc.symbols[index].side_size.x, module_count, x, y);
				}
				metadata_index += 3;
			}
		}
		//color palette
		for(jab_int32 i=2; i< std::min(enc.color_number, 64); i++)	//skip the first two colors in finder pattern
		{
			enc.symbols[index].matrix  [y*enc.symbols[index].side_size.x+x] = palette_index[master_palette_placement_index[0][i]%enc.color_number];
			enc.symbols[index].data_map[y*enc.symbols[index].side_size.x+x] = 0;
			module_count++;
			getNextMetadataModuleInMaster(enc.symbols[index].side_size.y, enc.symbols[index].side_size.x, module_count, x, y);

			enc.symbols[index].matrix  [y*enc.symbols[index].side_size.x+x] = palette_index[master_palette_placement_index[1][i]%enc.color_number];
			enc.symbols[index].data_map[y*enc.symbols[index].side_size.x+x] = 0;
			module_count++;
			getNextMetadataModuleInMaster(enc.symbols[index].side_size.y, enc.symbols[index].side_size.x, module_count, x, y);

			enc.symbols[index].matrix  [y*enc.symbols[index].side_size.x+x] = palette_index[master_palette_placement_index[2][i]%enc.color_number];
			enc.symbols[index].data_map[y*enc.symbols[index].side_size.x+x] = 0;
			module_count++;
			getNextMetadataModuleInMaster(enc.symbols[index].side_size.y, enc.symbols[index].side_size.x, module_count, x, y);

			enc.symbols[index].matrix  [y*enc.symbols[index].side_size.x+x] = palette_index[master_palette_placement_index[3][i]%enc.color_number];
			enc.symbols[index].data_map[y*enc.symbols[index].side_size.x+x] = 0;
			module_count++;
			getNextMetadataModuleInMaster(enc.symbols[index].side_size.y, enc.symbols[index].side_size.x, module_count, x, y);
		}
		//metadata PartII
		if(!isDefaultMode(enc))
		{
			while(metadata_index < enc.symbols[index].metadata.size())
			{
				color_index = 0;
				for(jab_int32 j=0; j<nb_of_bits_per_mod; j++)
				{
					if(metadata_index < enc.symbols[index].metadata.size())
					{
						color_index += ((jab_int32)enc.symbols[index].metadata[metadata_index]) << (nb_of_bits_per_mod-1-j);
						metadata_index++;
					}
					else
						break;
				}
				enc.symbols[index].matrix  [y*enc.symbols[index].side_size.x + x] = color_index;
				enc.symbols[index].data_map[y*enc.symbols[index].side_size.x + x] = 0;
				module_count++;
				getNextMetadataModuleInMaster(enc.symbols[index].side_size.y, enc.symbols[index].side_size.x, module_count, x, y);
			}
		}
    }
    else//place color palette in slave symbol
    {
    	//color palette
        jab_int32 width=enc.symbols[index].side_size.x;
        jab_int32 height=enc.symbols[index].side_size.y;
        for (jab_int32 i=2; i< std::min(enc.color_number, 64); i++)	//skip the first two colors in alignment pattern
        {
        	//left
			enc.symbols[index].matrix  [slave_palette_position[i-2].y*width + slave_palette_position[i-2].x] = palette_index[slave_palette_placement_index[i]%enc.color_number];
			enc.symbols[index].data_map[slave_palette_position[i-2].y*width + slave_palette_position[i-2].x] = 0;
			//top
			enc.symbols[index].matrix  [slave_palette_position[i-2].x*width + (width-1-slave_palette_position[i-2].y)] = palette_index[slave_palette_placement_index[i]%enc.color_number];
			enc.symbols[index].data_map[slave_palette_position[i-2].x*width + (width-1-slave_palette_position[i-2].y)] = 0;
			//right
			enc.symbols[index].matrix  [(height-1-slave_palette_position[i-2].y)*width + (width-1-slave_palette_position[i-2].x)] = palette_index[slave_palette_placement_index[i]%enc.color_number];
			enc.symbols[index].data_map[(height-1-slave_palette_position[i-2].y)*width + (width-1-slave_palette_position[i-2].x)] = 0;
			//bottom
			enc.symbols[index].matrix  [(height-1-slave_palette_position[i-2].x)*width + slave_palette_position[i-2].y] = palette_index[slave_palette_placement_index[i]%enc.color_number];
			enc.symbols[index].data_map[(height-1-slave_palette_position[i-2].x)*width + slave_palette_position[i-2].y] = 0;
        }
    }
    //Data placement
    jab_int32 written_mess_part=0;
    jab_int32 padding=0;
    for(jab_int32 start_i=0;start_i<enc.symbols[index].side_size.x;start_i++)
    {
        for(jab_int32 i=start_i;i<enc.symbols[index].side_size.x*enc.symbols[index].side_size.y;i=i+enc.symbols[index].side_size.x)
        {
            if (enc.symbols[index].data_map[i]!=0 && written_mess_part<ecc_encoded_data.size())
            {
                color_index=0;
                for(jab_int32 j=0;j<nb_of_bits_per_mod;j++)
                {
                    if(written_mess_part<ecc_encoded_data.size())
                        color_index+=((jab_int32)ecc_encoded_data[written_mess_part]) << (nb_of_bits_per_mod-1-j);//*pow(2,nb_of_bits_per_mod-1-j);
                    else
                    {
                        color_index+=padding << (nb_of_bits_per_mod-1-j);//*pow(2,nb_of_bits_per_mod-1-j);
                        if (padding==0)
                            padding=1;
                        else
                            padding=0;
                    }
                    written_mess_part++;
                }
                enc.symbols[index].matrix[i]=(jab_char)color_index;//i % enc.color_number;
            }
            else if(enc.symbols[index].data_map[i]!=0) //write padding bits
            {
                color_index=0;
                for(jab_int32 j=0;j<nb_of_bits_per_mod;j++)
                {
                    color_index+=padding << (nb_of_bits_per_mod-1-j);//*pow(2,nb_of_bits_per_mod-1-j);
                    if (padding==0)
                        padding=1;
                    else
                        padding=0;
                }
                enc.symbols[index].matrix[i]=(jab_char)color_index;//i % enc.color_number;

            }
        }
    }
	return JAB_SUCCESS;
}

/**
 * @brief Swap two symbols
 * @param enc the encode parameters
 * @param index1 the index number of the first symbol
 * @param index2 the index number of the second symbol
*/
constexpr static void swap_symbols(Encode& enc, jab_int32 index1, jab_int32 index2) noexcept {
	std::swap(enc.symbol_positions[index1],  enc.symbol_positions[index2]);
	std::swap(enc.symbol_versions[index1].x, enc.symbol_versions[index2].x);
	std::swap(enc.symbol_versions[index1].y, enc.symbol_versions[index2].y);
	std::swap(enc.symbol_ecc_levels[index1], enc.symbol_ecc_levels[index2]);
	//swap symbols
	std::swap(enc.symbols[index1], enc.symbols[index2]);
}

/**
 * @brief Assign docked symbols to their hosts
 * @param enc the encode parameters
 * @return JAB_SUCCESS | JAB_FAILURE
*/
constexpr static jab_boolean assignDockedSymbols(Encode& enc) noexcept {
	//initialize host and slaves
	for(jab_int32 i=0; i<enc.symbol_number; i++)
	{
		//initialize symbol host index
		enc.symbols[i].host = -1;
		//initialize symbol's slave index
		for(jab_int32 j=0; j<4; j++)
			enc.symbols[i].slaves[j] = 0;	//0:no slave
	}
	//assign docked symbols
	jab_int32 assigned_slave_index = 1;
    for(jab_int32 i=0; i<enc.symbol_number-1 && assigned_slave_index<enc.symbol_number; i++)
    {
    	for(jab_int32 j=0; j<4 && assigned_slave_index<enc.symbol_number; j++)
		{
			for(jab_int32 k=i+1; k<enc.symbol_number && assigned_slave_index<enc.symbol_number; k++)
			{
				if(enc.symbols[k].host == -1)
				{
					jab_int32 hpos = enc.symbol_positions[i];
					jab_int32 spos = enc.symbol_positions[k];
					jab_boolean slave_found = JAB_FAILURE;
					switch(j)
					{
					case 0:	//top
						if(jab_symbol_pos[hpos].x == jab_symbol_pos[spos].x && jab_symbol_pos[hpos].y - 1 == jab_symbol_pos[spos].y)
						{
							enc.symbols[i].slaves[0] = assigned_slave_index;
							enc.symbols[k].slaves[1] = -1;	//-1:host position
							slave_found = JAB_SUCCESS;
						}
						break;
					case 1:	//bottom
						if(jab_symbol_pos[hpos].x == jab_symbol_pos[spos].x && jab_symbol_pos[hpos].y + 1 == jab_symbol_pos[spos].y)
						{
							enc.symbols[i].slaves[1] = assigned_slave_index;
							enc.symbols[k].slaves[0] = -1;
							slave_found = JAB_SUCCESS;
						}
						break;
					case 2:	//left
						if(jab_symbol_pos[hpos].y == jab_symbol_pos[spos].y && jab_symbol_pos[hpos].x - 1 == jab_symbol_pos[spos].x)
						{
							enc.symbols[i].slaves[2] = assigned_slave_index;
							enc.symbols[k].slaves[3] = -1;
							slave_found = JAB_SUCCESS;
						}
						break;
					case 3://right
						if(jab_symbol_pos[hpos].y == jab_symbol_pos[spos].y && jab_symbol_pos[hpos].x + 1 == jab_symbol_pos[spos].x)
						{
							enc.symbols[i].slaves[3] = assigned_slave_index;
							enc.symbols[k].slaves[2] = -1;
							slave_found = JAB_SUCCESS;
						}
						break;
					}
					if(slave_found)
					{
						swap_symbols(enc, k, assigned_slave_index);
						enc.symbols[assigned_slave_index].host = i;
						assigned_slave_index++;
					}
				}
			}
		}
    }
    //check if there is undocked symbol
    for(jab_int32 i=1; i<enc.symbol_number; i++)
    {
		if(enc.symbols[i].host == -1)
		{
			JAB_REPORT_ERROR("Slave symbol at position ", enc.symbol_positions[i], " has no host");
			return JAB_FAILURE;
		}
    }
    return JAB_SUCCESS;
}

/**
 * @brief Calculate the code parameters according to the input symbols
 * @param enc the encode parameters
 * @return the code parameters
*/
constexpr static CodeParams getCodePara(Encode& enc) noexcept {
	CodeParams cp;
    //calculate the module size in pixel
    if(enc.master_symbol_width != 0 || enc.master_symbol_height != 0)
    {
        jab_int32 dimension_x = enc.master_symbol_width/enc.symbols[0].side_size.x;
        jab_int32 dimension_y = enc.master_symbol_height/enc.symbols[0].side_size.y;
        cp.dimension = dimension_x > dimension_y ? dimension_x : dimension_y;
		if (cp.dimension < 1) {
			cp.dimension = 1;
		}
    }
    else {
        cp.dimension = enc.module_size;
    }

    //find the coordinate range of symbols
    cp.min_x = 0;
    cp.min_y = 0;
    jab_int32 max_x=0, max_y=0;
    for(jab_int32 i=0; i<enc.symbol_number; i++)
    {
        //find the mininal x and y
		if (jab_symbol_pos[enc.symbol_positions[i]].x < cp.min_x) {
			cp.min_x = jab_symbol_pos[enc.symbol_positions[i]].x;
		}
		if (jab_symbol_pos[enc.symbol_positions[i]].y < cp.min_y) {
			cp.min_y = jab_symbol_pos[enc.symbol_positions[i]].y;
		}
        //find the maximal x and y
		if (jab_symbol_pos[enc.symbol_positions[i]].x > max_x) {
			max_x = jab_symbol_pos[enc.symbol_positions[i]].x;
		}
		if (jab_symbol_pos[enc.symbol_positions[i]].y > max_y) {
			max_y = jab_symbol_pos[enc.symbol_positions[i]].y;
		}
    }

    //calculate the code size
    cp.rows = max_y - cp.min_y + 1;
    cp.cols = max_x - cp.min_x + 1;
	cp.row_height = std::vector<jab_int32>(cp.rows);
    cp.col_width = std::vector<jab_int32>(cp.cols);
    cp.code_size.x = 0;
    cp.code_size.y = 0;
    jab_boolean flag = 0;
    for(jab_int32 x=cp.min_x; x<=max_x; x++)
    {
        flag = 0;
        for(jab_int32 i=0; i<enc.symbol_number; i++)
        {
            if(jab_symbol_pos[enc.symbol_positions[i]].x == x)
            {
                cp.col_width[x - cp.min_x] = enc.symbols[i].side_size.x;
                cp.code_size.x += cp.col_width[x - cp.min_x];
                flag = 1;
            }
            if(flag) break;
        }
    }
    for(jab_int32 y=cp.min_y; y<=max_y; y++)
    {
        flag = 0;
        for(jab_int32 i=0; i<enc.symbol_number; i++)
        {
            if(jab_symbol_pos[enc.symbol_positions[i]].y == y)
            {
                cp.row_height[y - cp.min_y] = enc.symbols[i].side_size.y;
                cp.code_size.y += cp.row_height[y - cp.min_y];
                flag = 1;
            }
            if(flag) break;
        }
    }
    return cp;
}

/**
 * @brief Create bitmap for the code
 * @param enc the encode parameters
 * @param cp the code parameters
 * @return JAB_SUCCESS | JAB_FAILURE
*/
constexpr static Bitmap createBitmap(Encode& enc, const CodeParams& cp) noexcept {
    //create bitmap
    jab_int32 width = cp.dimension * cp.code_size.x;
    jab_int32 height= cp.dimension * cp.code_size.y;
    jab_int32 bytes_per_pixel = BITMAP_BITS_PER_PIXEL / 8;
    jab_int32 bytes_per_row = width * bytes_per_pixel;
	Bitmap bp = Bitmap(
		width, height, BITMAP_BITS_PER_PIXEL, BITMAP_BITS_PER_CHANNEL, BITMAP_CHANNEL_COUNT,
		width*height*bytes_per_pixel, 0
	);

    //place symbols in bitmap
    for(jab_int32 k=0; k<enc.symbol_number; k++)
    {
        //calculate the starting coordinates of the symbol matrix
        jab_int32 startx = 0, starty = 0;
        jab_int32 col = jab_symbol_pos[enc.symbol_positions[k]].x - cp.min_x;
        jab_int32 row = jab_symbol_pos[enc.symbol_positions[k]].y - cp.min_y;
        for(jab_int32 c=0; c<col; c++)
            startx += cp.col_width[c];
        for(jab_int32 r=0; r<row; r++)
            starty += cp.row_height[r];

        //place symbol in the code
        jab_int32 symbol_width = enc.symbols[k].side_size.x;
        jab_int32 symbol_height= enc.symbols[k].side_size.y;
        for(jab_int32 x=startx; x<(startx+symbol_width); x++)
        {
            for(jab_int32 y=starty; y<(starty+symbol_height); y++)
            {
                //place one module in the bitmap
                jab_int32 p_index = enc.symbols[k].matrix[(y-starty)*symbol_width + (x-startx)];
                 for(jab_int32 i=y*cp.dimension; i<(y*cp.dimension+cp.dimension); i++)
                {
                    for(jab_int32 j=x*cp.dimension; j<(x*cp.dimension+cp.dimension); j++)
                    {
                        bp.pixel[i*bytes_per_row + j*bytes_per_pixel]     = enc.palette[p_index * 3];	//R
                        bp.pixel[i*bytes_per_row + j*bytes_per_pixel + 1] = enc.palette[p_index * 3 + 1];//B
                        bp.pixel[i*bytes_per_row + j*bytes_per_pixel + 2] = enc.palette[p_index * 3 + 2];//G
                        bp.pixel[i*bytes_per_row + j*bytes_per_pixel + 3] = 255; 							//A
                    }
                }
            }
        }
    }
    return bp;
}


/**
 * @brief Checks if the docked symbol sizes are valid
 * @param enc the encode parameters
 * @return JAB_SUCCESS | JAB_FAILURE
*/
constexpr static jab_boolean checkDockedSymbolSize(const Encode& enc) noexcept {
	for(jab_int32 i=0; i<enc.symbol_number; i++)
	{
		for(jab_int32 j=0; j<4; j++)
		{
			jab_int32 slave_index = enc.symbols[i].slaves[j];
			if(slave_index > 0)
			{
				jab_int32 hpos = enc.symbol_positions[i];
				jab_int32 spos = enc.symbol_positions[slave_index];
				jab_int32 x_diff = jab_symbol_pos[hpos].x - jab_symbol_pos[spos].x;
				jab_int32 y_diff = jab_symbol_pos[hpos].y - jab_symbol_pos[spos].y;

				if(x_diff == 0 && enc.symbol_versions[i].x != enc.symbol_versions[slave_index].x)
				{
					JAB_REPORT_ERROR("Slave symbol at position ", spos, "has different side version in X direction as its host symbol at position ", hpos);
					return JAB_FAILURE;
				}
				if(y_diff == 0 && enc.symbol_versions[i].y != enc.symbol_versions[slave_index].y)
				{
					JAB_REPORT_ERROR("Slave symbol at position ", spos, "has different side version in Y direction as its host symbol at position", hpos);
					return JAB_FAILURE;
				}
			}
		}
	}
	return JAB_SUCCESS;
}

/**
 * @brief Set the minimal master symbol version
 * @param enc the encode parameters
 * @param encoded_data the encoded message
 * @return JAB_SUCCESS | JAB_FAILURE
 */
static jab_boolean setMasterSymbolVersion(Encode& enc, std::vector<jab_char>& encoded_data) noexcept {
    //calculate required number of data modules depending on data_length
    jab_int32 net_data_length = encoded_data.size();
    jab_int32 payload_length = net_data_length + 5;  //plus S and flag bit
    if(enc.symbol_ecc_levels[0] == 0) enc.symbol_ecc_levels[0] = DEFAULT_ECC_LEVEL;
    enc.symbols[0].wcwr[0] = ecclevel2wcwr[enc.symbol_ecc_levels[0]][0];
	enc.symbols[0].wcwr[1] = ecclevel2wcwr[enc.symbol_ecc_levels[0]][1];

	//determine the minimum square symbol to fit data
	jab_int32 capacity, net_capacity;
	jab_boolean found_flag = JAB_FAILURE;
	for (jab_int32 i=1; i<=32; i++)
	{
		enc.symbol_versions[0].x = i;
		enc.symbol_versions[0].y = i;
		capacity = getSymbolCapacity(enc, 0);
		net_capacity = (capacity/enc.symbols[0].wcwr[1])*enc.symbols[0].wcwr[1] - (capacity/enc.symbols[0].wcwr[1])*enc.symbols[0].wcwr[0];
		if(net_capacity >= payload_length)
		{
			found_flag = JAB_SUCCESS;
			break;
		}
	}
	if(!found_flag)
	{
		jab_int32 level = -1;
		for (jab_int32 j=(jab_int32)enc.symbol_ecc_levels[0]-1; j>0; j--)
		{
			net_capacity = (capacity/ecclevel2wcwr[j][1])*ecclevel2wcwr[j][1] - (capacity/ecclevel2wcwr[j][1])*ecclevel2wcwr[j][0];
			if(net_capacity >= payload_length)
				level = j;
		}
		if(level > 0)
		{
			JAB_REPORT_ERROR("Message does not fit into one symbol with the given ECC level. Please use an ECC level lower than ", level);
			return JAB_FAILURE;
		}
		else
		{
			reportError("Message does not fit into one symbol. Use more symbols.");
			return JAB_FAILURE;
		}
	}
	//update symbol side size
    enc.symbols[0].side_size.x = VERSION2SIZE(enc.symbol_versions[0].x);
	enc.symbols[0].side_size.y = VERSION2SIZE(enc.symbol_versions[0].y);

    return JAB_SUCCESS;
}

/**
 * @brief Add variable E to slave symbol metadata the data payload for each symbol
 * @param slave the slave symbol
 * @return JAB_SUCCESS | JAB_FAILURE
*/
static jab_boolean addE2SlaveMetadata(Symbol* slave) noexcept {
	//copy old metadata to new metadata
	jab_int32 old_metadata_length = slave->metadata.size();
	jab_int32 new_metadata_length = old_metadata_length + 6;
	std::vector<jab_char> old_metadata = std::move(slave->metadata);
	slave->metadata = std::vector<jab_char>(new_metadata_length);

	memcpy(slave->metadata.data(), old_metadata.data(), old_metadata_length);

	//update SE = 1
	slave->metadata[1] = 1;
	//set variable E
	jab_int32 E1 = slave->wcwr[0] - 3;
	jab_int32 E2 = slave->wcwr[1] - 4;
	convert_dec_to_bin(E1, slave->metadata.data(), old_metadata_length, 3);
	convert_dec_to_bin(E2, slave->metadata.data(), old_metadata_length + 3, 3);
	return JAB_SUCCESS;
}

/**
 * @brief Update slave metadata E in its host data stream
 * @param enc the encode parameters
 * @param host_index the host symbol index
 * @param slave_index the slave symbol index
*/
static void updateSlaveMetadataE(Encode& enc, jab_int32 host_index, jab_int32 slave_index) noexcept {
	Symbol& host = enc.symbols[host_index];
	const Symbol& slave = enc.symbols[slave_index];

	jab_int32 offset = host.data.size() - 1;
	//find the start flag of metadata
	while (host.data[offset] == 0)
	{
		offset--;
	}
	//skip the flag bit
	offset--;
	//skip host metadata S
	if (host_index == 0)
		offset -= 4;
	else
		offset -= 3;
	//skip other slave symbol's metadata
	for (jab_int32 i = 0; i < 4; i++)
	{
		if (host.slaves[i] == slave_index)
			break;
		else if (host.slaves[i] <= 0)
			continue;
		else
			offset -= enc.symbols[host.slaves[i]].metadata.size();
	}
	//skip SS, SE and possibly V
	if (slave.metadata[0] == 1) {
		offset -= 7;
	}
	else {
		offset -= 2;
	}

	//update E
	std::array<jab_char, 6> E = {};
	jab_int32 E1 = slave.wcwr[0] - 3;
	jab_int32 E2 = slave.wcwr[1] - 4;
	convert_dec_to_bin(E1, E.data(), 0, 3);
	convert_dec_to_bin(E2, E.data(), 3, 3);
	for(jab_int32 i=0; i<6; i++)
	{
		host.data[offset--] = E[i];
	}
}

/**
 * @brief Set the data payload for each symbol
 * @param enc the encode parameters
 * @param encoded_data the encoded message
 * @return JAB_SUCCESS | JAB_FAILURE
*/
constexpr static jab_boolean fitDataIntoSymbols(Encode& enc, std::vector<jab_char>& encoded_data) noexcept {
	//calculate the net capacity of each symbol and the total net capacity
	std::vector<jab_int32> capacity = std::vector<jab_int32>(enc.symbol_number);
	std::vector<jab_int32> net_capacity = std::vector<jab_int32>(enc.symbol_number);
	jab_int32 total_net_capacity = 0;
	for(jab_int32 i = 0; i < enc.symbol_number; i++) {
		capacity[i] = getSymbolCapacity(enc, i);
		enc.symbols[i].wcwr[0] = ecclevel2wcwr[enc.symbol_ecc_levels[i]][0];
		enc.symbols[i].wcwr[1] = ecclevel2wcwr[enc.symbol_ecc_levels[i]][1];
		net_capacity[i] = (capacity[i]/enc.symbols[i].wcwr[1])*enc.symbols[i].wcwr[1] - (capacity[i]/enc.symbols[i].wcwr[1])*enc.symbols[i].wcwr[0];
		total_net_capacity += net_capacity[i];
	}
	//assign data into each symbol
	jab_int32 assigned_data_length = 0;
	for(jab_int32 i=0; i<enc.symbol_number; i++)
	{
		//divide data proportionally
		jab_int32 s_data_length;
		if(i == enc.symbol_number - 1)
		{
			s_data_length = encoded_data.size() - assigned_data_length;
		}
		else
		{
			jab_float prop = (jab_float)net_capacity[i] / (jab_float)total_net_capacity;
			s_data_length = (jab_int32)(prop * encoded_data.size());
		}
		jab_int32 s_payload_length = s_data_length;

		//add flag bit
		s_payload_length++;
		//add host metadata S length (master metadata Part III or slave metadata Part III)
		if(i == 0)	s_payload_length += 4;
		else		s_payload_length += 3;
		//add slave metadata length
		for(jab_int32 j=0; j<4; j++)
		{
			if(enc.symbols[i].slaves[j] > 0)
			{
				s_payload_length += enc.symbols[enc.symbols[i].slaves[j]].metadata.size();
			}
		}

		//check if the full payload exceeds net capacity
		if(s_payload_length > net_capacity[i])
		{
			reportError("Message does not fit into the specified code. Use higher symbol version.");
			return JAB_FAILURE;
		}

		//add metadata E for slave symbols if free capacity available
		jab_int32 j = 0;
		while(net_capacity[i] - s_payload_length >= 6 && j < 4)
		{
			if(enc.symbols[i].slaves[j] > 0)
			{
				if(enc.symbols[enc.symbols[i].slaves[j]].metadata[1] == 0) //check SE
				{
					if (!addE2SlaveMetadata(&enc.symbols[enc.symbols[i].slaves[j]])) {
						return JAB_FAILURE;
					}
					s_payload_length += 6;	//add E length
				}
			}
			j++;
		}

		//get optimal code rate
		jab_int32 pn_length = s_payload_length;
		if(i == 0)
		{
			if(!isDefaultMode(enc))
			{
				getOptimalECC(capacity[i], s_payload_length, enc.symbols[i].wcwr);
				pn_length = (capacity[i]/enc.symbols[i].wcwr[1])*enc.symbols[i].wcwr[1] - (capacity[i]/enc.symbols[i].wcwr[1])*enc.symbols[i].wcwr[0];
			}
			else
				pn_length = net_capacity[i];
		}
		else
		{
			if(enc.symbols[i].metadata[1] == 1)	//SE = 1
			{
				getOptimalECC(capacity[i], pn_length, enc.symbols[i].wcwr);
				pn_length = (capacity[i]/enc.symbols[i].wcwr[1])*enc.symbols[i].wcwr[1] - (capacity[i]/enc.symbols[i].wcwr[1])*enc.symbols[i].wcwr[0];
				updateSlaveMetadataE(enc, enc.symbols[i].host, i);
			}
			else
				pn_length = net_capacity[i];
		}

		//start to set full payload
        enc.symbols[i].data = std::vector<jab_char>(pn_length);
		//set data
		memcpy(enc.symbols[i].data.data(), &encoded_data[assigned_data_length], s_data_length);
		assigned_data_length += s_data_length;
		//set flag bit
		jab_int32 set_pos = s_payload_length - 1;
		enc.symbols[i].data[set_pos--] = 1;
		//set host metadata S
		for(jab_int32 k=0; k<4; k++)
		{
			if(enc.symbols[i].slaves[k] > 0)
			{
				enc.symbols[i].data[set_pos--] = 1;
			}
			else if(enc.symbols[i].slaves[k] == 0)
			{
				enc.symbols[i].data[set_pos--] = 0;
			}
		}
		//set slave metadata
		for(jab_int32 k=0; k<4; k++)
		{
			if(enc.symbols[i].slaves[k] > 0)
			{
				for(jab_int32 m=0; m<enc.symbols[enc.symbols[i].slaves[k]].metadata.size(); m++)
				{
					enc.symbols[i].data[set_pos--] = enc.symbols[enc.symbols[i].slaves[k]].metadata[m];
				}
			}
		}
	}
	return JAB_SUCCESS;
}

/**
 * @brief Initialize symbols
 * @param enc the encode parameters
 * @return JAB_SUCCESS | JAB_FAILURE
*/
constexpr static jab_boolean InitSymbols(Encode& enc) noexcept {
	//check all information for multi-symbol code are valid
	if(enc.symbol_number > 1)
	{
		for(jab_int32 i=0; i<enc.symbol_number; i++)
		{
			if(enc.symbol_versions[i].x < 1 || enc.symbol_versions[i].x > 32 || enc.symbol_versions[i].y < 1 || enc.symbol_versions[i].y > 32)
			{
				JAB_REPORT_ERROR("Incorrect symbol version for symbol ", i);
				return JAB_FAILURE;
			}
			if(enc.symbol_positions[i] < 0 || enc.symbol_positions[i] > MAX_SYMBOL_NUMBER)
			{
				JAB_REPORT_ERROR("Incorrect symbol position for symbol ", i);
				return JAB_FAILURE;
			}
		}
	}
	//move the master symbol to the first
    if(enc.symbol_number > 1 && enc.symbol_positions[0] != 0)
	{
		for(jab_int32 i=0; i<enc.symbol_number; i++)
		{
			if(enc.symbol_positions[i] == 0)
			{
				std::swap(enc.symbol_positions[i], 	enc.symbol_positions[0]);
				std::swap(enc.symbol_versions[i].x, 	enc.symbol_versions[0].x);
				std::swap(enc.symbol_versions[i].y, 	enc.symbol_versions[0].y);
				std::swap(enc.symbol_ecc_levels[i], 	enc.symbol_ecc_levels[0]);
				break;
			}
		}
	}
    //if no master symbol exists in multi-symbol code
    if(enc.symbol_number > 1 && enc.symbol_positions[0] != 0)
    {
		reportError("Master symbol missing");
		return JAB_FAILURE;
    }
    //if only one symbol but its position is not 0 - set to zero. Everything else makes no sense.
	if (enc.symbol_number == 1 && enc.symbol_positions[0] != 0) {
		enc.symbol_positions[0] = 0;
	}
    //check if a symbol position is used twice
	for(jab_int32 i=0; i<enc.symbol_number-1; i++)
	{
		for(jab_int32 j=i+1; j<enc.symbol_number; j++)
		{
			if(enc.symbol_positions[i] == enc.symbol_positions[j])
			{
				reportError("Duplicate symbol position");
				return JAB_FAILURE;
			}
		}
	}
	//assign docked symbols to their hosts
	if (!assignDockedSymbols(enc)) {
		return JAB_FAILURE;
	}
    //check if the docked symbol size matches the docked side of its host
	if (!checkDockedSymbolSize(enc)) {
		return JAB_FAILURE;
	}
	//set symbol index and symbol side size
    for(jab_int32 i=0; i<enc.symbol_number; i++)
    {
    	//set symbol index
		enc.symbols[i].index = i;
        //set symbol side size
        enc.symbols[i].side_size.x = VERSION2SIZE(enc.symbol_versions[i].x);
        enc.symbols[i].side_size.y = VERSION2SIZE(enc.symbol_versions[i].y);
    }
    return JAB_SUCCESS;
}

/**
 * @brief Set metadata for slave symbols
 * @param enc the encode parameters
 * @return JAB_SUCCESS | JAB_FAILURE
*/
constexpr static jab_boolean setSlaveMetadata(Encode& enc) noexcept {
	//set slave metadata variables
	for(jab_int32 i=1; i<enc.symbol_number; i++)
	{
		jab_int32 SS, SE, V, E1=0, E2=0;
		jab_int32 metadata_length = 2; //Part I length
		//SS and V
		if(enc.symbol_versions[i].x != enc.symbol_versions[enc.symbols[i].host].x)
		{
		  	SS = 1;
			V = enc.symbol_versions[i].x - 1;
			metadata_length += 5;
		}
		else if(enc.symbol_versions[i].y != enc.symbol_versions[enc.symbols[i].host].y)
		{
			SS = 1;
			V = enc.symbol_versions[i].y - 1;
			metadata_length += 5;
		}
		else
		{
			SS = 0;
		}
		//SE and E
		if(enc.symbol_ecc_levels[i] == 0 || enc.symbol_ecc_levels[i] == enc.symbol_ecc_levels[enc.symbols[i].host])
		{
			SE = 0;
		}
		else
		{
			SE = 1;
			E1 = ecclevel2wcwr[enc.symbol_ecc_levels[i]][0] - 3;
			E2 = ecclevel2wcwr[enc.symbol_ecc_levels[i]][1] - 4;
			metadata_length += 6;
		}
		//write slave metadata
		enc.symbols[i].metadata = std::vector<jab_char>(metadata_length);
		//Part I
		enc.symbols[i].metadata[0] = SS;
		enc.symbols[i].metadata[1] = SE;
		//Part II
		if(SS == 1)
		{
			convert_dec_to_bin(V, enc.symbols[i].metadata.data(), 2, 5);
		}
		if(SE == 1)
		{
			jab_int32 start_pos = (SS == 1) ? 7 : 2;
			convert_dec_to_bin(E1, enc.symbols[i].metadata.data(), start_pos, 3);
			convert_dec_to_bin(E2, enc.symbols[i].metadata.data(), start_pos + 3, 3);
		}
	}
	return JAB_SUCCESS;
}

/**
 * @brief Generate JABCode
 * @param enc the encode parameters
 * @param data the input data
 * @return 0:success | 1: out of memory | 2:no input data | 3:incorrect symbol version or position | 4: input data too long
*/
std::optional<Bitmap> generateJABCode(Encode& enc, const std::vector<jab_char>& data, jab_int32& error) noexcept {
    //Check data
    if(data.size() == 0)
    {
        reportError("No input data specified!");
		error = 2;
        return std::nullopt;
    }

    //initialize symbols and set metadata in symbols
	if (!InitSymbols(enc)) {
		error = 3;
		return std::nullopt;
	}

    //get the optimal encoded length and encoding sequence
    jab_int32 encoded_length;
    std::optional<std::vector<jab_int32>> encode_seq = analyzeInputData(data, encoded_length);
    if(!encode_seq.has_value())
	{
		reportError("Analyzing input data failed");
		error = 1;
		return std::nullopt;
    }
	//encode data using optimal encoding modes
    std::optional<std::vector<jab_char>> encoded_dataOp = encodeData(data, encoded_length, encode_seq.value());
    if(!encoded_dataOp.has_value())
    {
		reportError("Error encoding data");
		error = 1;
		return std::nullopt;
    }

	std::vector<jab_char> encoded_data = std::move(encoded_dataOp.value());
    //set master symbol version if not given
    if(enc.symbol_number == 1 && (enc.symbol_versions[0].x == 0 || enc.symbol_versions[0].y == 0))
    {
        if(!setMasterSymbolVersion(enc, encoded_data))
        {
			error = 4;
            return std::nullopt;
        }
    }
	//set metadata for slave symbols
	if(!setSlaveMetadata(enc))
	{
		error = 1;
		return std::nullopt;
	}
	//assign encoded data into symbols
	if(!fitDataIntoSymbols(enc, encoded_data))
	{
		error = 4;
		return std::nullopt;
	}
	//set master metadata
	if(!isDefaultMode(enc))
	{
		if(!encodeMasterMetadata(enc))
		{
			reportError("Encoding master symbol metadata failed");
			error = 1;
            return std::nullopt;
		}
	}

    //encode each symbol in turn
    for(jab_int32 i=0; i<enc.symbol_number; i++)
    {
        //error correction for data
        std::vector<jab_char> ecc_encoded_data = encodeLDPC(enc.symbols[i].data.data(), enc.symbols[i].wcwr, enc.symbols[i].data.size());

        //interleave
        interleaveData(ecc_encoded_data);
        //create Matrix
        jab_boolean cm_flag = createMatrix(enc, i, ecc_encoded_data);
        if(!cm_flag)
        {
			JAB_REPORT_ERROR("Creating matrix for symbol ", i, " failed");
			error = 1;
			return std::nullopt;
		}
    }

    //mask all symbols in the code
    CodeParams cp = getCodePara(enc);
    if(isDefaultMode(enc))	//default mode
	{
		maskSymbols(enc, DEFAULT_MASKING_REFERENCE);
	}
	else
	{
		jab_int32 mask_reference = maskCode(enc, cp);
		if(mask_reference < 0)
		{
			error = 1;
			return std::nullopt;
		}
		if(mask_reference != DEFAULT_MASKING_REFERENCE)
		{
			//re-encode PartII of master symbol metadata
			updateMasterMetadataPartII(enc, mask_reference);
			//update the masking reference in master symbol metadata
			placeMasterMetadataPartII(enc);
		}
	}

    //create the code bitmap
    return createBitmap(enc, cp);
}
