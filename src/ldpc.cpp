
#include "ldpc.h"

#include "reportError.h"

#include <cmath>
#include "detector.h"
#include "pseudo_random.h"
#include <memory>

/**
 * @brief Create matrix A for message data
 * @param wc the number of '1's in a column
 * @param wr the number of '1's in a row
 * @param capacity the number of columns of the matrix
 * @return the matrix A | NULL if failed (out of memory)
*/
static std::unique_ptr<jab_int32[]> createMatrixA(jab_int32 wc, jab_int32 wr, jab_int32 capacity) noexcept {
    jab_int32 nb_pcb = wr < 4 ? capacity / 2 : capacity / wr * wc;

    jab_int32 effwidth = std::ceil(capacity/32.0f) * 32;
    jab_int32 offset = std::ceil(capacity/32.0f);
    //create a matrix with '0' entries
    std::unique_ptr<jab_int32[]> matrixA = std::unique_ptr<jab_int32[]>(new jab_int32[std::ceil(capacity / 32.0f) * nb_pcb]{});

    //Fill the first set with consecutive ones in each row
    for (jab_int32 y = 0; y < capacity / wr; y++) {
        for (jab_int32 x = 0; x < wr; x++) {
            matrixA[(y * (effwidth + wr) + x) / 32] |= 1 << (31 - ((y * (effwidth + wr) + x) % 32));
        }
    }

    jab_int32* permutation = new jab_int32[capacity];
    for (jab_int32 i = 0; i < capacity; i++) {
        permutation[i] = i;
    }

    //Permutate the columns and fill the remaining matrix
    //generate matrixA by following Gallagers algorithm
    LCG lcg = LCG(LPDC_MESSAGE_SEED);
    for (jab_int32 i = 1; i < wc; i++) {
        jab_int32 off_index = i * (capacity / wr);
        for (jab_int32 j = 0; j < capacity; j++) {
            jab_int32 pos = (jab_int32)((jab_float)lcg.temper() / (jab_float)UINT32_MAX * (capacity - j));
            for (jab_int32 k = 0; k < capacity / wr; k++) {
                matrixA[(off_index + k) * offset + j / 32] |= ((matrixA[(permutation[pos] / 32 + k * offset)] >> (31 - permutation[pos] % 32)) & 1) << (31 - j % 32);
            }
            std::swap(permutation[capacity - 1 - j], permutation[pos]);
        }
    }
    delete[] permutation;
    return matrixA;
}

#include <memory>
static enum class EncodeMode {
    Encode,
    Decode
};
/**
 * @brief Gauss Jordan elimination algorithm
 * @param matrixA the matrix
 * @param wc the number of '1's in a column
 * @param wr the number of '1's in a row
 * @param capacity the number of columns of the matrix
 * @param matrix_rank the rank of the matrix
 * @param encode specifies if function is called by the encoder or decoder (true is encode)
 * @return 0: success | 1: fatal error (out of memory)
*/
template<EncodeMode mode>
constexpr static void GaussJordan(std::unique_ptr<jab_int32[]>& matrixA, jab_int32 wc, jab_int32 wr, jab_int32 capacity, jab_int32& matrix_rank) noexcept {
    jab_int32 loop=0;
    jab_int32 nb_pcb = wr < 4 ? capacity / 2 : capacity/wr*wc;

    jab_int32 offset = std::ceil(capacity/(jab_float)32);

    // It is significantly faster with raw pointers compared to std::vector or std::unique_ptr
    jab_int32* column_arrangement = new jab_int32[capacity]{};
    jab_boolean* processed_column = new jab_boolean[capacity]{};
    jab_int32* zero_lines_nb = new jab_int32[nb_pcb]{};
    jab_int32* swap_col = new jab_int32[2 * capacity];

    jab_int32* matrixH = new jab_int32[offset * nb_pcb]{};
    std::memcpy(matrixH, matrixA.get(), offset * nb_pcb * sizeof(jab_int32));

    jab_int32 zero_lines=0;

    for (jab_int32 i=0; i<nb_pcb; i++)
    {
        jab_int32 pivot_column=capacity+1;
        for (jab_int32 j=0; j<capacity; j++)
        {
            if((matrixH[(offset*32*i+j)/32] >> (31-(offset*32*i+j)%32)) & 1)
            {
                pivot_column=j;
                break;
            }
        }
        if(pivot_column < capacity)
        {
            processed_column[pivot_column]=1;
            column_arrangement[pivot_column]=i;
            if (pivot_column>=nb_pcb)
            {
                swap_col[2*loop]=pivot_column;
                loop++;
            }

            jab_int32 off_index=pivot_column/32;
            jab_int32 off_index1=pivot_column%32;
            for (jab_int32 j=0; j<nb_pcb; j++)
            {
                if (j != i && ((matrixH[off_index+j*offset] >> (31-off_index1)) & 1))
                {
                    //subtract pivot row GF(2)
                    for (jab_int32 k = 0; k < offset; k++) {
                        matrixH[k + offset * j] ^= matrixH[k + offset * i];
                    }
                }
            }
        }
        else //zero line
        {
            zero_lines_nb[zero_lines]=i;
            zero_lines++;
        }
    }

    matrix_rank = nb_pcb-zero_lines;
    jab_int32 loop2=0;
    for(jab_int32 i = matrix_rank; i < nb_pcb; i++)
    {
        if(column_arrangement[i] > 0)
        {
            for (jab_int32 j=0;j < nb_pcb;j++)
            {
                if (processed_column[j] == 0)
                {
                    column_arrangement[j]=column_arrangement[i];
                    column_arrangement[i]=0;
                    processed_column[j]=1;
                    processed_column[i]=0;
                    swap_col[2*loop]=i;
                    swap_col[2*loop+1]=j;
                    column_arrangement[i]=j;
                    loop++;
                    loop2++;
                    break;
                }
            }
        }
    }

    jab_int32 loop1=0;
    for (jab_int32 kl=0;kl< nb_pcb;kl++)
    {
        if((loop1 < loop - loop2) && processed_column[kl] == 0)
        {
            column_arrangement[kl]=column_arrangement[swap_col[2*loop1]];
            processed_column[kl]=1;
            swap_col[2*loop1+1]=kl;
            loop1++;
        }
    }

    loop1=0;
    for (jab_int32 kl=0;kl< nb_pcb;kl++)
    {
        if(processed_column[kl]==0)
        {
            column_arrangement[kl]=zero_lines_nb[loop1];
            loop1++;
        }
    }
    //rearrange matrixH if encoder and store it in matrixA
    //rearrange matrixA if decoder
    if constexpr (mode == EncodeMode::Encode) {
        for(jab_int32 i=0;i< nb_pcb;i++)
            std::memcpy(matrixA.get() + i * offset, matrixH + column_arrangement[i] * offset, offset * sizeof(jab_int32));

        //swap columns
        jab_int32 tmp=0;
#define invertedLoops1 false // Originally it is false
#if !invertedLoops1
        for(jab_int32 i=0;i<loop;i++)
        {
            for (jab_int32 j=0;j<nb_pcb;j++)
            {
                tmp ^= (-((matrixA[swap_col[2*i]/32+j*offset] >> (31-swap_col[2*i]%32)) & 1) ^ tmp) & (1 << 0);
                matrixA[swap_col[2*i]/32+j*offset]   ^= (-((matrixA[swap_col[2*i+1]/32+j*offset] >> (31-swap_col[2*i+1]%32)) & 1) ^ matrixA[swap_col[2*i]/32+j*offset]) & (1 << (31-swap_col[2*i]%32));
                matrixA[swap_col[2*i+1]/32+offset*j] ^= (-((tmp >> 0) & 1) ^ matrixA[swap_col[2*i+1]/32+offset*j]) & (1 << (31-swap_col[2*i+1]%32));
            }
        }
#else
        for (jab_int32 y = 0; y < nb_pcb; y++) {
            for (jab_int32 x = 0; x < loop; x++) {
                tmp ^= (-((matrixA[swap_col[2 * x] / 32 + y * offset] >> (31 - swap_col[2 * x] % 32)) & 1) ^ tmp) & (1 << 0);
                matrixA[swap_col[2 * x] / 32 + y * offset] ^= (-((matrixA[swap_col[2 * x + 1] / 32 + y * offset] >> (31 - swap_col[2 * x + 1] % 32)) & 1) ^ matrixA[swap_col[2 * x] / 32 + y * offset]) & (1 << (31 - swap_col[2 * x] % 32));
                matrixA[swap_col[2 * x + 1] / 32 + offset * y] ^= (-((tmp >> 0) & 1) ^ matrixA[swap_col[2 * x + 1] / 32 + offset * y]) & (1 << (31 - swap_col[2 * x + 1] % 32));
            }
        }
#endif
    }
    else {
        for (jab_int32 i = 0; i < nb_pcb; i++) {
            std::memcpy(matrixH + i * offset, matrixA.get() + column_arrangement[i] * offset, offset * sizeof(jab_int32));
        }

        //swap columns
        jab_int32 tmp=0;
#if !invertedLoops1
        for (jab_int32 i = 0; i < loop; i++)
        {
            for (jab_int32 j = 0; j < nb_pcb; j++)
            {
                tmp ^= (-((matrixH[swap_col[2*i]/32+j*offset] >> (31-swap_col[2*i]%32)) & 1) ^ tmp) & (1 << 0);
                matrixH[swap_col[2*i]/32+j*offset]   ^= (-((matrixH[swap_col[2*i+1]/32+j*offset] >> (31-swap_col[2*i+1]%32)) & 1) ^ matrixH[swap_col[2*i]/32+j*offset]) & (1 << (31-swap_col[2*i]%32));
                matrixH[swap_col[2*i+1]/32+offset*j] ^= (-((tmp >> 0) & 1) ^ matrixH[swap_col[2*i+1]/32+offset*j]) & (1 << (31-swap_col[2*i+1]%32));
            }
        }
#else
        for (jab_int32 y = 0; y < nb_pcb; y++) {
            for (jab_int32 x = 0; x < loop; x++) {
                tmp ^= (-((matrixH[swap_col[2 * x] / 32 + y * offset] >> (31 - swap_col[2 * x] % 32)) & 1) ^ tmp) & (1 << 0);
                matrixH[swap_col[2 * x] / 32 + y * offset] ^= (-((matrixH[swap_col[2 * x + 1] / 32 + y * offset] >> (31 - swap_col[2 * x + 1] % 32)) & 1) ^ matrixH[swap_col[2 * x] / 32 + y * offset]) & (1 << (31 - swap_col[2 * x] % 32));
                matrixH[swap_col[2 * x + 1] / 32 + offset * y] ^= (-((tmp >> 0) & 1) ^ matrixH[swap_col[2 * x + 1] / 32 + offset * y]) & (1 << (31 - swap_col[2 * x + 1] % 32));
            }
        }
#endif
        std::memcpy(matrixA.get(), matrixH, offset* nb_pcb * sizeof(jab_int32));
    }
    delete[] column_arrangement;
    delete[] processed_column;
    delete[] zero_lines_nb;
    delete[] swap_col;
    delete[] matrixH;
}

/**
 * @brief Create the error correction matrix for the metadata
 * @param wc the number of '1's in a column
 * @param capacity the number of columns of the matrix
 * @return the error correction matrix | NULL if failed
*/
static std::unique_ptr<jab_int32[]> createMetadataMatrixA(jab_int32 wc, jab_int32 capacity) noexcept {
    jab_int32 nb_pcb = capacity/2;
    jab_int32 offset = std::ceil(capacity/(jab_float)32);
    //create a matrix with '0' entries


    jab_int32* permutation = new jab_int32[capacity];
    for (jab_int32 i = 0; i < capacity; i++) {
        permutation[i] = i;
    }

    LCG lcg = LCG(LPDC_METADATA_SEED);
    jab_int32 nb_once=capacity*nb_pcb/(jab_float)wc+3;
    nb_once=nb_once/nb_pcb;

    std::unique_ptr<jab_int32[]> matrixA = std::unique_ptr<jab_int32[]>(new jab_int32[offset * nb_pcb]{});
    //Fill matrix randomly
    for (jab_int32 i=0;i<nb_pcb;i++)
    {
        for (jab_int32 j=0; j< nb_once; j++)
        {
            jab_int32 pos = (jab_int32)( (jab_float)lcg.temper() / (jab_float)UINT32_MAX * (capacity-j) );
            matrixA[i*offset+permutation[pos]/32] |= 1 << (31-permutation[pos]%32);
            std::swap(permutation[capacity - 1 - j], permutation[pos]);
        }
    }
    delete[] permutation;
    return matrixA;
}

/**
 * @brief Create the generator matrix to encode messages
 * @param matrixA the error correction matrix
 * @param capacity the number of columns of the matrix
 * @param Pn the number of net message bits
 * @return the generator matrix | NULL if failed (out of memory)
*/
static std::unique_ptr<jab_int32[]> createGeneratorMatrix(const std::unique_ptr<jab_int32[]>& matrixA, jab_int32 capacity, jab_int32 Pn) noexcept {
    jab_int32 effwidth=std::ceil(Pn/(jab_float)32)*32;
    jab_int32 offset=std::ceil(Pn/(jab_float)32);
    jab_int32 offset_cap=std::ceil(capacity/(jab_float)32);
    //create G [I C]
    //remember matrixA is now A = [I CT], now use it and create G=[CT
                //                                                 I ]

    std::unique_ptr<jab_int32[]> G = std::unique_ptr<jab_int32[]>(new jab_int32[offset * capacity]{});
    //fill identity matrix
    for (jab_int32 i = 0; i < Pn; i++) {
        G[(capacity - Pn + i) * offset + i / 32] |= 1 << (31 - i % 32);
    }

    //copy CT matrix from A to G
    jab_int32 matrix_index=capacity-Pn;
    jab_int32 loop=0;

    for (jab_int32 i=0; i<(capacity-Pn)*effwidth; i++)
    {
        if(matrix_index >= capacity)
        {
            loop++;
            matrix_index=capacity-Pn;
        }
        if(i % effwidth < Pn)
        {
            G[i/32] ^= (-((matrixA[matrix_index/32+offset_cap*loop] >> (31-matrix_index%32)) & 1) ^ G[i/32]) & (1 << (31-i%32));
            matrix_index++;
        }
    }
    return G;
}

/**
 * @brief LDPC encoding
 * @param data the data to be encoded
 * @param coderate_params the two code rate parameter wc and wr indicating how many '1' in a column (Wc) and how many '1' in a row of the parity check matrix
 * @return the encoded data | NULL if failed
*/
std::vector<jab_char> encodeLDPC(jab_char* data, jab_int32* coderate_params, const jab_int32 size) noexcept {
    jab_int32 matrix_rank=0;
    jab_int32 wc, wr, Pg, Pn;       //number of '1' in column //number of '1' in row //gross message length //number of parity check symbols //calculate required parameters
    wc = coderate_params[0];
    wr = coderate_params[1];
    Pn = size;
    if(wr > 0) {
        Pg = std::ceil((Pn*wr)/(jab_float)(wr-wc));
        Pg = wr * (std::ceil(Pg / (jab_float)wr));
    }
    else {
        Pg = Pn * 2;
    }


    //in order to speed up the ldpc encoding, sub blocks are established
    jab_int32 nb_sub_blocks=0;
    for(jab_int32 i = 1;i<10000;i++)
    {
        if(Pg / i < 2700)
        {
            nb_sub_blocks=i;
            break;
        }
    }
    jab_int32 Pg_sub_block=0;
    jab_int32 Pn_sub_block=0;
    if(wr > 0)
    {
        Pg_sub_block=((Pg / nb_sub_blocks) / wr) * wr;
        Pn_sub_block=Pg_sub_block * (wr-wc) / wr;
    }
    else
    {
        Pg_sub_block=Pg;
        Pn_sub_block=Pn;
    }
    jab_int32 encoding_iterations=nb_sub_blocks=Pg / Pg_sub_block;//nb_sub_blocks;
    if (Pn_sub_block * nb_sub_blocks < Pn) {
        encoding_iterations--;
    }

    std::unique_ptr<jab_int32[]> matrixA = wr > 0 ? createMatrixA(wc, wr, Pg_sub_block) :
                                              createMetadataMatrixA(wc, Pg_sub_block);

    GaussJordan<EncodeMode::Encode>(matrixA, wc, wr, Pg_sub_block, matrix_rank);

    //Generator Matrix
    std::unique_ptr<jab_int32[]> G = createGeneratorMatrix(matrixA, Pg_sub_block, Pg_sub_block - matrix_rank);

    std::vector<jab_char> ecc_encoded_data = std::vector<jab_char>(Pg);

    jab_int32 temp,loop;
    jab_int32 offset=std::ceil((Pg_sub_block - matrix_rank)/(jab_float)32);
    //G * message = ecc_encoded_Data
    for(jab_int32 iter=0; iter < encoding_iterations; iter++)
    {
        for (jab_int32 i=0;i<Pg_sub_block;i++)
        {
            temp=0;
            loop=0;
            jab_int32 offset_index=offset*i;
            for (jab_int32 j=iter*Pn_sub_block; j < (iter+1)*Pn_sub_block; j++)
            {
                temp ^= (((G[offset_index + loop/32] >> (31-loop%32)) & 1) & ((data[j] >> 0) & 1)) << 0;
                loop++;
            }
            ecc_encoded_data[i+iter*Pg_sub_block]=(jab_char) ((temp >> 0) & 1);
        }
    }
    if(encoding_iterations != nb_sub_blocks)
    {
        jab_int32 start=encoding_iterations*Pn_sub_block;
        jab_int32 last_index=encoding_iterations*Pg_sub_block;
        matrix_rank=0;
        Pg_sub_block=Pg - encoding_iterations * Pg_sub_block;
        Pn_sub_block=Pg_sub_block * (wr-wc) / wr;

        std::unique_ptr<jab_int32[]> matrixA2 = createMatrixA(wc, wr, Pg_sub_block);

        GaussJordan<EncodeMode::Encode>(matrixA2, wc, wr, Pg_sub_block, matrix_rank);
        std::unique_ptr<jab_int32[]> G = createGeneratorMatrix(matrixA2, Pg_sub_block, Pg_sub_block - matrix_rank);

        offset = std::ceil((Pg_sub_block - matrix_rank)/(jab_float)32);
        for (jab_int32 i=0;i<Pg_sub_block;i++)
        {
            temp=0;
            loop=0;
            jab_int32 offset_index=offset*i;
            for (jab_int32 j=start; j < size; j++)
            {
                temp ^= (((G[offset_index + loop/32] >> (31-loop%32)) & 1) & ((data[j] >> 0) & 1)) << 0;
                loop++;
            }
            ecc_encoded_data[i+last_index]=(jab_char) ((temp >> 0) & 1);
        }
    }
    return ecc_encoded_data;
}

/**
 * @brief Iterative hard decision error correction decoder
 * @param data the received data
 * @param matrix the parity check matrix
 * @param length the encoded data length
 * @param height the number of check bits
 * @param max_iter the maximal number of iterations
 * @param is_correct indicating if decodedMessage function could correct all errors
 * @param start_pos indicating the position to start reading in data array
 * @return 1: error correction succeeded | 0: fatal error (out of memory)
*/
static jab_int32 decodeMessage(jab_byte* data, const jab_int32* matrix, jab_int32 length, jab_int32 height, jab_int32 max_iter, jab_boolean& is_correct, jab_int32 start_pos) noexcept {
    jab_int32* max_val = new jab_int32[length]{};
    jab_int32* equal_max = new jab_int32[length]{};
    jab_int32* prev_index = new jab_int32[length]{};

    is_correct=(jab_boolean)1;
    jab_int32 check=0;
    jab_int32 counter=0, prev_count=0;
    jab_int32 max=0;
    jab_int32 offset=std::ceil(length/(jab_float)32);

    for (jab_int32 kl=0;kl<max_iter;kl++)
    {
        max=0;
        for(jab_int32 j=0;j<height;j++)
        {
            check=0;
            for (jab_int32 i=0;i<length;i++)
            {
                if(((matrix[j*offset+i/32] >> (31-i%32)) & 1) & ((data[start_pos+i] >> 0) & 1))
                    check+=1;
            }
            check=check%2;
            if(check)
            {
                for(jab_int32 k=0;k<length;k++)
                {
                    if(((matrix[j*offset+k/32] >> (31-k%32)) & 1))
                        max_val[k]++;
                }
            }
        }
        //find maximal values in max_val
        jab_boolean is_used=0;
        for (jab_int32 j=0;j<length;j++)
        {
            is_used=(jab_boolean)0;
            for(jab_int32 i=0;i< prev_count;i++)
            {
                if(prev_index[i]==j)
                    is_used=(jab_boolean)1;
            }
            if(max_val[j]>=max && !is_used)
            {
                if(max_val[j]!=max)
                    counter=0;
                max=max_val[j];
                equal_max[counter]=j;
                counter++;
            }
            max_val[j]=0;
        }
        //flip bits
        if(max>0)
        {
            is_correct=(jab_boolean) 0;
            if(length < 36)
            {
                jab_int32 rand_tmp=(jab_int32)(rand()/(jab_float)UINT32_MAX * counter);
                prev_index[0]=start_pos+equal_max[rand_tmp];
                data[start_pos+equal_max[rand_tmp]]=(data[start_pos+equal_max[rand_tmp]]+1)%2;
            }
            else
            {
                for(jab_int32 j=0; j< counter;j++)
                {
                    prev_index[j]=start_pos+equal_max[j];
                    data[start_pos+equal_max[j]]=(data[start_pos+equal_max[j]]+1)%2;
                }
            }
            prev_count=counter;
            counter=0;
        }
        else
            is_correct=(jab_boolean) 1;

        if(is_correct == 0 && kl+1 < max_iter)
            is_correct=(jab_boolean)1;
        else
            break;
    }
    delete[] max_val;
    delete[] equal_max;
    delete[] prev_index;
    return 1;
}

/**
 * @brief LDPC decoding to perform hard decision
 * @param data the encoded data
 * @param length the encoded data length
 * @param wc the number of '1's in a column
 * @param wr the number of '1's in a row
 * @return the decoded data length | 0: fatal error (out of memory)
*/
jab_int32 decodeLDPChd(jab_byte* data, jab_int32 length, jab_int32 wc, jab_int32 wr) noexcept {
    jab_int32 matrix_rank=0;
    jab_int32 max_iter=25;
    jab_int32 Pn, Pg, decoded_data_len = 0;
    if(wr > 3)
    {
        Pg = wr * (length / wr);
        Pn = Pg * (wr - wc) / wr;                //number of source symbols
    }
    else
    {
        Pg=length;
        Pn=length/2;
        wc=2;
        if(Pn>36)
            wc=3;
    }
    decoded_data_len=Pn;

    //in order to speed up the ldpc encoding, sub blocks are established
    jab_int32 nb_sub_blocks=0;
    jab_int32 i = 1;
    for(;i<10000;i++)
    {
        if(Pg / i < 2700)
        {
            nb_sub_blocks=i;
            break;
        }
    }
    jab_int32 Pg_sub_block=0;
    jab_int32 Pn_sub_block=0;
    if(wr > 3)
    {
        Pg_sub_block=((Pg / nb_sub_blocks) / wr) * wr;
        Pn_sub_block=Pg_sub_block * (wr-wc) / wr;
    }
    else
    {
        Pg_sub_block=Pg;
        Pn_sub_block=Pn;
    }
    jab_int32 decoding_iterations=nb_sub_blocks=Pg / Pg_sub_block;//nb_sub_blocks;
    if(Pn_sub_block * nb_sub_blocks < Pn)
        decoding_iterations--;

    //parity check matrix
    std::unique_ptr<jab_int32[]> matrixA = wr > 0 ? createMatrixA(wc, wr, Pg_sub_block) : createMetadataMatrixA(wc, Pg_sub_block);

    GaussJordan<EncodeMode::Decode>(matrixA, wc, wr, Pg_sub_block, matrix_rank);

    jab_int32 old_Pg_sub=Pg_sub_block;
    jab_int32 old_Pn_sub=Pn_sub_block;
    for (jab_int32 iter = 0; iter < nb_sub_blocks; iter++)
    {
        if(decoding_iterations != nb_sub_blocks && iter == decoding_iterations)
        {
            matrix_rank=0;
            Pg_sub_block=Pg - decoding_iterations * Pg_sub_block;
            Pn_sub_block=Pg_sub_block * (wr-wc) / wr;

            std::unique_ptr<jab_int32[]> matrixA1 = createMatrixA(wc, wr, Pg_sub_block);

            GaussJordan<EncodeMode::Decode>(matrixA1, wc, wr, Pg_sub_block, matrix_rank);

            //ldpc decoding
            //first check syndrom
            jab_boolean is_correct=1;
            jab_int32 offset=std::ceil(Pg_sub_block/(jab_float)32);
            for (jab_int32 i=0;i< matrix_rank; i++)
            {
                jab_int32 temp=0;
                for (jab_int32 j=0;j<Pg_sub_block;j++)
                    temp ^= (((matrixA1[i*offset+j/32] >> (31-j%32)) & 1) & ((data[iter*old_Pg_sub+j] >> 0) & 1)) << 0; //
                if (temp)
                {
                    is_correct=(jab_boolean) 0;//message not correct
                    break;
                }
            }

            if(is_correct==0) {
                jab_int32 start_pos=iter*old_Pg_sub;
                jab_int32 success = decodeMessage(data, matrixA1.get(), Pg_sub_block, matrix_rank, max_iter, is_correct, start_pos);
                if(success == 0)
                {
                    reportError("LDPC decoder error.");
                    return 0;
                }
            }
            if(is_correct==0)
            {
                jab_boolean is_correct=1;
                for (jab_int32 i=0;i< matrix_rank; i++)
                {
                    jab_int32 temp=0;
                    for (jab_int32 j = 0; j < Pg_sub_block; j++) {
                        temp ^= (((matrixA1[i * offset + j / 32] >> (31 - j % 32)) & 1) & ((data[iter * old_Pg_sub + j] >> 0) & 1)) << 0;
                    }
                    if (temp)
                    {
                        is_correct=(jab_boolean) 0;//message not correct
                        break;
                    }
                }
                if(is_correct==0)
                {
                    reportError("Too many errors in message. LDPC decoding failed.");
                    return 0;
                }
            }
        }
        else {
            //ldpc decoding
            //first check syndrom
            jab_boolean is_correct=1;
            jab_int32 offset=std::ceil(Pg_sub_block/(jab_float)32);
            for (jab_int32 i=0;i< matrix_rank; i++)
            {
                jab_int32 temp=0;
                for (jab_int32 j=0;j<Pg_sub_block;j++)
                    temp ^= (((matrixA[i*offset+j/32] >> (31-j%32)) & 1) & ((data[iter*old_Pg_sub+j] >> 0) & 1)) << 0;
                if (temp)
                {
                    is_correct=(jab_boolean) 0;//message not correct
                    break;
                }
            }

            if(is_correct==0)
            {
                jab_int32 start_pos=iter*old_Pg_sub;
                jab_int32 success = decodeMessage(data, matrixA.get(), Pg_sub_block, matrix_rank, max_iter, is_correct, start_pos);
                if(success == 0)
                {
                    reportError("LDPC decoder error.");
                    return 0;
                }
                is_correct=1;
                for (jab_int32 i=0;i< matrix_rank; i++)
                {
                    jab_int32 temp=0;
                    for (jab_int32 j=0;j<Pg_sub_block;j++)
                        temp ^= (((matrixA[i*offset+j/32] >> (31-j%32)) & 1) & ((data[iter*old_Pg_sub+j] >> 0) & 1)) << 0;
                    if (temp)
                    {
                        is_correct=(jab_boolean)0;//message not correct
                        break;
                    }
                }
                if(is_correct==0)
                {
                    reportError("Too many errors in message. LDPC decoding failed.");
                    return 0;
                }
            }
        }
        jab_int32 loop=0;
        for (jab_int32 i=iter*old_Pg_sub;i < iter * old_Pg_sub + Pn_sub_block; i++)
        {
            data[iter*old_Pn_sub+loop]=data[i+ matrix_rank];
            loop++;
        }
    }
    return decoded_data_len;
}

/**
 * @brief LDPC Iterative belief propagation decoding algorithm for binary codes
 * @param enc the received reliability value for each bit
 * @param matrix the decoding matrixreliability value for each bit
 * @param length the encoded data length
 * @param checkbits the rank of the matrix
 * @param height the number of check bits
 * @param max_iter the maximal number of iterations
 * @param is_correct indicating if decodedMessage function could correct all errors
 * @param start_pos indicating the position to start reading in enc array
 * @param dec is the tentative decision after each decoding iteration
 * @return 1: error correction succeded | 0: decoding failed
*/
static jab_int32 decodeMessageBP(jab_float* enc, const jab_int32* matrix, jab_int32 length, jab_int32 checkbits, jab_int32 height, jab_int32 max_iter, jab_boolean *is_correct, jab_int32 start_pos, jab_byte* dec) noexcept {
    jab_double* lambda = new jab_double[length];
    jab_double* old_nu_row = new jab_double[length];
    jab_double* nu = new jab_double[length * height]{};
    jab_int32* index = new jab_int32[length];

    jab_int32 offset=std::ceil(length/(jab_float)32);
    jab_double product=1.0;

    //set last bits
    for (jab_int32 i=length-1;i >= length-(height-checkbits);i--)
    {
        enc[start_pos+i]=1.0;
        dec[start_pos+i]=0;
    }

    jab_double meansum=0.0;
    for (jab_int32 i=0;i<length;i++)
        meansum+=enc[start_pos+i];

    //calc variance
    meansum/=length;
    jab_double var=0.0;
    for (jab_int32 i=0;i<length;i++)
        var+=(enc[start_pos+i]-meansum)*(enc[start_pos+i]-meansum);
    var/=(length-1);

    //initialize lambda
    for (jab_int32 i=0;i<length;i++)
    {
        if(dec[start_pos+i])
            enc[start_pos+i]=-enc[start_pos+i];
        lambda[i]=(jab_double)2.0*enc[start_pos+i]/var;
    }

    //check node update
    jab_int32 count;
    for (jab_int32 kl=0;kl<max_iter;kl++)
    {
        for(jab_int32 j=0;j<height;j++)
        {
            product=1.0;
            count=0;
            for (jab_int32 i=0;i<length;i++)
            {
                if((matrix[j*offset+i/32] >> (31-i%32)) & 1)
                {
                    if (kl==0)
                        product*=std::tanh(lambda[i]*0.5);
                    else
                        product*=std::tanh(nu[j*length+i]*0.5);
                    index[count]=i;
                    ++count;
                }
            }
            //update nu
            jab_double num=0.0, denum=0.0;
            for (jab_int32 i=0;i<count;i++)
            {
                if(((matrix[j*offset+index[i]/32] >> (31-index[i]%32)) & 1) && tanh(nu[j*length+index[i]]*0.5) != 0.0 && kl>0)// && tanh(-(lambda[index[i]]-old_nu_row[index[i]])/2) != 0.0)
                {
                    num     = 1 + product / std::tanh(nu[j*length+index[i]]*0.5);
                    denum   = 1 - product / std::tanh(nu[j*length+index[i]]*0.5);
                }
                else if(((matrix[j*offset+index[i]/32] >> (31-index[i]%32)) & 1) && tanh(lambda[index[i]]*0.5) != 0.0 && kl==0)// && tanh(-(lambda[index[i]]-old_nu_row[index[i]])/2) != 0.0)
                {
                    num     = 1 + product / std::tanh(lambda[index[i]]*0.5);
                    denum   = 1 - product / std::tanh(lambda[index[i]]*0.5);
                }
                else
                {
                    num     = 1 + product;
                    denum   = 1 - product;
                }
                if (num == 0.0)
                    nu[j*length+index[i]]=-1;
                else if(denum == 0.0)
                    nu[j*length+index[i]]= 1;
                else
                    nu[j*length+index[i]]= std::log(num / denum);
            }
        }
        //update lambda
        jab_double sum = 0.0;
        for (jab_int32 i=0;i<length;i++) {
            for(jab_int32 k=0;k<height;k++)
            {
                sum+=nu[k*length+i];
                old_nu_row[k]=nu[k*length+i];
            }
            for(jab_int32 k=0;k<height;k++)
            {
                if((matrix[k*offset+i/32] >> (31-i%32)) & 1)
                    nu[k*length+i]=lambda[i]+(sum-old_nu_row[k]);
            }
            lambda[i]=2.0*enc[start_pos+i]/var+sum;
            if (lambda[i] < 0) {
                dec[start_pos + i] = 1;
            }
            else {
                dec[start_pos + i] = 0;
            }
            sum = 0.0;
        }
        //check matrix times dec
        *is_correct=(jab_boolean) 1;
        for (jab_int32 i=0;i< height; i++)
        {
            jab_int32 temp=0;
            for (jab_int32 j=0;j<length;j++)
                temp ^= (((matrix[i*offset+j/32] >> (31-j%32)) & 1) & ((dec[start_pos+j] >> 0) & 1)) << 0;
            if (temp)
            {
                *is_correct=(jab_boolean) 0;
                break;//message not correct
            }
        }
        if(!*is_correct && kl<max_iter-1)
            *is_correct=(jab_boolean) 1;
        else
            break;
    }

    delete[] lambda;
    delete[] old_nu_row;
    delete[] nu;
    delete[] index;
    return 1;
}

/**
 * @brief LDPC decoding to perform soft decision
 * @param enc the probability value for each bit position
 * @param length the encoded data length
 * @param wc the number of '1's in each column
 * @param wr the number of '1's in each row
 * @param dec the decoded data
 * @return the decoded data length | 0: decoding error
*/
jab_int32 decodeLDPC(jab_float* enc, jab_int32 length, jab_int32 wc, jab_int32 wr, jab_byte* dec) noexcept {
    jab_int32 matrix_rank=0;
    jab_int32 max_iter=25;
    jab_int32 Pn, Pg, decoded_data_len = 0;
    if(wr > 3)
    {
        Pg = wr * (length / wr);
        Pn = Pg * (wr - wc) / wr; //number of source symbols
    }
    else
    {
        Pg=length;
        Pn=length/2;
        wc=2;
        if (Pn > 36) {
            wc = 3;
        }
    }
    decoded_data_len=Pn;

    //in order to speed up the ldpc encoding, sub blocks are established
    jab_int32 nb_sub_blocks=0;
    for(jab_int32 i=1;i<10000;i++)
    {
        if(Pg / i < 2700)
        {
            nb_sub_blocks=i;
            break;
        }
    }
    jab_int32 Pg_sub_block=0;
    jab_int32 Pn_sub_block=0;
    if(wr > 3)
    {
        Pg_sub_block=((Pg / nb_sub_blocks) / wr) * wr;
        Pn_sub_block=Pg_sub_block * (wr-wc) / wr;
    }
    else
    {
        Pg_sub_block=Pg;
        Pn_sub_block=Pn;
    }
    jab_int32 decoding_iterations=nb_sub_blocks=Pg / Pg_sub_block;//nb_sub_blocks;
    if (Pn_sub_block * nb_sub_blocks < Pn) {
        decoding_iterations--;
    }

    //parity check matrix
    std::unique_ptr<jab_int32[]> matrixA = wr > 0 ? createMatrixA(wc, wr, Pg_sub_block) : createMetadataMatrixA(wc, Pg_sub_block);
    GaussJordan<EncodeMode::Decode>(matrixA, wc, wr, Pg_sub_block, matrix_rank);

    jab_int32 old_Pg_sub=Pg_sub_block;
    jab_int32 old_Pn_sub=Pn_sub_block;
    for (jab_int32 iter = 0; iter < nb_sub_blocks; iter++)
    {
        if(decoding_iterations != nb_sub_blocks && iter == decoding_iterations)
        {
            matrix_rank=0;
            Pg_sub_block=Pg - decoding_iterations * Pg_sub_block;
            Pn_sub_block=Pg_sub_block * (wr-wc) / wr;

            std::unique_ptr<jab_int32[]> matrixA1 = createMatrixA(wc, wr, Pg_sub_block);
            GaussJordan<EncodeMode::Decode>(matrixA1, wc, wr, Pg_sub_block, matrix_rank);

            //ldpc decoding
            //first check syndrom
            jab_boolean is_correct=1;
            jab_int32 offset=std::ceil(Pg_sub_block/(jab_float)32);
            for (jab_int32 i=0;i< matrix_rank; i++)
            {
                jab_int32 temp=0;
                for (jab_int32 j=0;j<Pg_sub_block;j++)
                    temp ^= (((matrixA1[i*offset+j/32] >> (31-j%32)) & 1) & ((dec[iter*old_Pg_sub+j] >> 0) & 1)) << 0; //
                if (temp)
                {
                    is_correct=(jab_boolean) 0; //message not correct
                    break;
                }
            }

            if(is_correct==0)
            {
                jab_int32 start_pos=iter*old_Pg_sub;
                jab_int32 success = decodeMessageBP(enc, matrixA1.get(), Pg_sub_block, matrix_rank, wr<4 ? Pg_sub_block / 2 : Pg_sub_block / wr * wc, max_iter, & is_correct, start_pos, dec);
                if(success == 0)
                {
                    reportError("LDPC decoder error.");
                    return 0;
                }
            }
            if(is_correct==0)
            {
                jab_boolean is_correct=1;
                for (jab_int32 i=0;i< matrix_rank; i++)
                {
                    jab_int32 temp=0;
                    for (jab_int32 j=0;j<Pg_sub_block;j++)
                        temp ^= (((matrixA1[i*offset+j/32] >> (31-j%32)) & 1) & ((dec[iter*old_Pg_sub+j] >> 0) & 1)) << 0;
                    if (temp)
                    {
                        is_correct=(jab_boolean) 0;//message not correct
                        break;
                    }
                }
                if(is_correct==0)
                {
                    reportError("Too many errors in message. LDPC decoding failed.");
                    return 0;
                }
            }
        }
        else
        {
            //ldpc decoding
            //first check syndrom
            jab_boolean is_correct=1;
            jab_int32 offset=std::ceil(Pg_sub_block/(jab_float)32);
            for (jab_int32 i=0;i< matrix_rank; i++)
            {
                jab_int32 temp=0;
                for (jab_int32 j=0;j<Pg_sub_block;j++)
                    temp ^= (((matrixA[i*offset+j/32] >> (31-j%32)) & 1) & ((dec[iter*old_Pg_sub+j] >> 0) & 1)) << 0;
                if (temp)
                {
                    is_correct=(jab_boolean) 0;//message not correct
                    break;
                }
            }

            if(is_correct==0)
            {
                jab_int32 start_pos=iter*old_Pg_sub;
                jab_int32 success = decodeMessageBP(enc, matrixA.get(), Pg_sub_block, matrix_rank, wr<4 ? Pg_sub_block / 2 : Pg_sub_block / wr * wc, max_iter, & is_correct, start_pos, dec);
                if(success == 0)
                {
                    reportError("LDPC decoder error.");
                    return 0;
                }
                is_correct=1;
                for (jab_int32 i=0;i< matrix_rank; i++)
                {
                    jab_int32 temp=0;
                    for (jab_int32 j=0;j<Pg_sub_block;j++)
                        temp ^= (((matrixA[i*offset+j/32] >> (31-j%32)) & 1) & ((dec[iter*old_Pg_sub+j] >> 0) & 1)) << 0;
                    if (temp)
                    {
                        is_correct=(jab_boolean)0;//message not correct
                        break;
                    }
                }
                if(is_correct==0)
                {
                    reportError("Too many errors in message. LDPC decoding failed.");
                    return 0;
                }
            }
        }

        jab_int32 loop=0;
        for (jab_int32 i=iter*old_Pg_sub;i < iter * old_Pg_sub + Pn_sub_block; i++) {
            dec[iter*old_Pn_sub+loop]=dec[i+ matrix_rank];
            loop++;
        }
    }
    return decoded_data_len;
}
