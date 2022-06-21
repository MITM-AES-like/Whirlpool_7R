#include "rndfill.h"
#include "aes.h"
#include "defines.h"

#define Round_N (7)
#define DoM_InBytes (4)
#define DoF_BL_InBytes (4)
#define DoF_RD_InBytes (1)
#define GUESS_BL_InBytes (0)
#define GUESS_RD_InBytes (3)
#define GUESS_BR_InBytes (0)
#define DoM_InBits (DoM_InBytes * SB_n)
#define DoF_BL_InBits (DoF_BL_InBytes * SB_n)
#define DoF_RD_InBits (DoF_RD_InBytes * SB_n)
#define GUESS_BL_InBits (GUESS_BL_InBytes * SB_n)
#define GUESS_RD_InBits (GUESS_RD_InBytes * SB_n)
#define GUESS_BR_InBits (GUESS_BR_InBytes * SB_n)

// For full-state match, CONST_SAMPLE_N should be (1 << ((STATE_InBytes - (DoF_BL_InBytes + DoF_RD_InBytes)) * SB_n))
#define CONST_SAMPLE_N (1ULL << (PARTIAL_InBits - (DoFF + DoFB) + 8)) // 1ULL << (PARTIAL_InBits - (DoFF + DoFB)) is sufficient, + 8 to leave some redundance repetitions

using namespace std;

string fn[omp_nb_threads];
ofstream fout[omp_nb_threads];

// Fix the main diagonal of AK3 and k3, such that
// the main diagonal of SB4_C be iS[0], such that after SB and SR, the first column of AK4 is 0
// After adding k4, the first column of SB5 equals the first column of k4
// Such that, the anti-main diagonal of MC5 equals that of KMC4.
// After adding MC5 with KMC4, the Red cells cancelled. All cells in SB6 are Blue.
const state_t k3_fixed_C = _mm_set_epi8(State_to_Parameter_Order(
    iS[0],    _ ,    _ ,   _ ,
       _ , iS[0],    _ ,   _ ,
       _ ,    _ , iS[0],   _ ,
       _ ,    _ ,    _ , iS[0]));
// Fix the cells in AK3[0, 3, 5, 10] such that one can compute the 
// initial values of Blue cells in MC3, such that all the values have constant impact on the
// main diagonal of SB4
const state_t AK3_C = _mm_set_epi8(State_to_Parameter_Order(
     0x00,    _ ,    _ , 0x00,
       _ ,  0x00,    _ , 0x00,
       _ ,    _ ,  0x00, 0x00,
     0x00,    _ ,    _ , 0x00));
const state_t AK3_C_MASK = _mm_set_epi8(State_to_Parameter_Order(
    0xFF,   _ ,   _ , 0xFF,
      _ , 0xFF,   _ , 0xFF,
      _ ,   _ , 0xFF, 0xFF,
    0xFF,   _ ,   _ , 0xFF));

const state_t SB3_C = _mm_set_epi8(State_to_Parameter_Order(
      _ ,   _ ,   _ , 0x00,
    0x00,   _ ,   _ ,   _ ,
      _ , 0x00,   _ ,   _ ,
      _ ,   _ , 0x00,   _ ));
const state_t SB3_C_MASK = _mm_set_epi8(State_to_Parameter_Order(
      _ ,   _ ,   _ , 0xFF,
    0xFF,   _ ,   _ ,   _ ,
      _ , 0xFF,   _ ,   _ ,
      _ ,   _ , 0xFF,   _ ));
const state_t SB3_RD_MASK = SB3_C_MASK;
const state_t SB3_BL_MASK = _mm_xor_si128(SB3_C_MASK, _mm_set1_epi8(0xFF));

const state_t MC5_C = _mm_set_epi8(State_to_Parameter_Order(
    0x00,   _ ,   _ ,   _ ,
      _ ,   _ ,   _ , 0x00,
      _ ,   _ , 0x00,   _ ,
      _ , 0x00,   _ ,   _ ));
const state_t MC5_RD_MASK = _mm_set_epi8(State_to_Parameter_Order(
    0xFF,   _ ,   _ ,   _ ,
      _ ,   _ ,   _ , 0xFF,
      _ ,   _ , 0xFF,   _ ,
      _ , 0xFF,   _ ,   _ ));
const state_t MC5_BL_MASK = _mm_xor_si128(MC5_RD_MASK, _mm_set1_epi8(0xFF));

const state_t MC2_C = _mm_set_epi8(State_to_Parameter_Order(
      _ ,   _ ,   _ , 0x00,
      _ ,   _ , 0x00,   _ ,
      _ , 0x00,   _ ,   _ ,
    0x00,   _ ,   _ ,   _ )); // Constants in #MC2_C
const state_t MC2_C_MASK = _mm_set_epi8(State_to_Parameter_Order(
      _ ,   _ ,   _ , 0xFF,
      _ ,   _ , 0xFF,   _ ,
      _ , 0xFF,   _ ,   _ ,
    0xFF,   _ ,   _ ,   _ ));

const state_t AK6_MATCH_MASK = _mm_set_epi8(State_to_Parameter_Order(
    0xFF,   _ ,   _ ,   _ ,
      _ , 0xFF,   _ ,   _ ,
      _ ,   _ , 0xFF,   _ ,
      _ ,   _ ,   _ , 0xFF));

const state_t AK0_RD_MASK = _mm_set_epi8(State_to_Parameter_Order(
      _ ,   _ ,   _ ,   _ ,
    0xFF,   _ ,   _ ,   _ ,
      _ ,   _ ,   _ ,   _ ,
      _ ,   _ ,   _ ,   _ ));
const state_t AK0_GUESS_MASK = _mm_set_epi8(State_to_Parameter_Order(
    0xFF,   _ ,   _ ,   _ ,
      _ ,   _ ,   _ ,   _ ,
    0xFF,   _ ,   _ ,   _ ,
    0xFF,   _ ,   _ ,   _ ));

const cell_t * AK3_C_pt = (const cell_t *)&AK3_C; // pointer to Constants in #AK3_C
const cell_t * MC2_C_pt = (const cell_t *)&MC2_C; // pointer to Constants in #MC2_C
const cell_t * SB3_C_pt = (const cell_t *)&SB3_C; // pointer to Constants in #SB3_C
const cell_t * MC5_C_pt = (const cell_t *)&MC5_C; // pointer to Constants in #MC5_C


map<cell_t, vector<array<cell_t, 4> > > AK3_RD_MC5_RD;
array<int, 1ULL<<24> counter = {0};
/**************************************************************************************
 * Precompute lookup table AK3_C_MC5_R,
 * the index is the values of the 3-byte impacts on #AK3[0, 5, 10],
 * the value is the Red cells (neutral bytes for Backward) in state #MC5
 * Because of the equivalent attack configuration, 
 * this precomputation is not used in the attack
 * Thus, this is only used for double check the degree of freedom of backward
 **************************************************************************************/
void testRD1()
{
    cout << "Start precompute the Red cells (neutral bytes for Backward) in state #MC5." << endl;
    uint64_t complexity = 0ULL;
    state_t MC5_RD  = _mm_setzero_si128();
    state_t KMC4_RD = _mm_setzero_si128();
    state_t SB5_RD  = _mm_setzero_si128();
    state_t k4_RD   = _mm_setzero_si128();
    state_t AK5_RD  = _mm_setzero_si128();
    state_t MC4_RD  = _mm_setzero_si128();
    state_t KMC3_RD = _mm_setzero_si128();
    state_t SB4_RD  = _mm_setzero_si128();
    state_t k3_RD   = _mm_setzero_si128();
    state_t AK3_RD  = _mm_setzero_si128();
    state_t KMC2_RD = _mm_setzero_si128();

    cell_t * MC5_RD_pt  = (cell_t * )&MC5_RD ;
    cell_t * KMC4_RD_pt = (cell_t * )&KMC4_RD;
    cell_t * SB5_RD_pt  = (cell_t * )&SB5_RD ;
    cell_t * k4_RD_pt   = (cell_t * )&k4_RD  ;
    cell_t * AK5_RD_pt  = (cell_t * )&AK5_RD ;
    cell_t * MC4_RD_pt  = (cell_t * )&MC4_RD ;
    cell_t * KMC3_RD_pt = (cell_t * )&KMC3_RD;
    cell_t * SB4_RD_pt  = (cell_t * )&SB4_RD ;
    cell_t * k3_RD_pt   = (cell_t * )&k3_RD  ;
    cell_t * AK3_RD_pt  = (cell_t * )&AK3_RD ;
    cell_t * KMC2_RD_pt = (cell_t * )&KMC2_RD;

    complexity = 0ULL;
    array<cell_t, 3> AK3_RD_C_0_5_10 = {0};
    array<cell_t, 4> MC5_RD_0_7_10_13 = {0};

    for (uint64_t i = 0ULL; i < (1ULL << 32); i++)
    {
        MC5_RD  = _mm_setzero_si128();
        KMC4_RD = _mm_setzero_si128();
        SB5_RD  = _mm_setzero_si128();
        k4_RD   = _mm_setzero_si128();
        AK5_RD  = _mm_setzero_si128();
        MC4_RD  = _mm_setzero_si128();
        KMC3_RD = _mm_setzero_si128();
        SB4_RD  = _mm_setzero_si128();
        k3_RD   = _mm_setzero_si128();
        AK3_RD  = _mm_setzero_si128();
        KMC2_RD = _mm_setzero_si128();

        MC5_RD_pt[ 0] = (cell_t)( i        & 0xffU);
        MC5_RD_pt[ 7] = (cell_t)((i >>  8) & 0xffU);
        MC5_RD_pt[10] = (cell_t)((i >> 16) & 0xffU);
        MC5_RD_pt[13] = (cell_t)((i >> 24) & 0xffU);

        KMC4_RD_pt[ 0] = (cell_t)( i        & 0xffU);
        KMC4_RD_pt[ 7] = (cell_t)((i >>  8) & 0xffU);
        KMC4_RD_pt[10] = (cell_t)((i >> 16) & 0xffU);
        KMC4_RD_pt[13] = (cell_t)((i >> 24) & 0xffU);

        SB5_RD = inv_SB_SR(MC5_RD);
        k4_RD = inv_SB_SR(KMC4_RD);

        SB4_RD = dec_one_round(SB5_RD, k4_RD);
        k3_RD = dec_one_round(k4_RD, RC[4]);

        AK3_RD = add_round_key(SB4_RD, k3_RD);

        counter[(((uint32_t)(AK3_RD_pt[ 0])) | ((uint32_t)(AK3_RD_pt[ 5]) << 8) | ((uint32_t)(AK3_RD_pt[10]) << 16))]++;

        if ((AK3_RD_pt[ 0] == AK3_C_pt[ 0]) && 
            (AK3_RD_pt[ 5] == AK3_C_pt[ 5]) && 
            (AK3_RD_pt[10] == AK3_C_pt[10]))
        {
            MC5_RD_0_7_10_13[0] = MC5_RD_pt[ 0];
            MC5_RD_0_7_10_13[1] = MC5_RD_pt[ 7];
            MC5_RD_0_7_10_13[2] = MC5_RD_pt[10];
            MC5_RD_0_7_10_13[3] = MC5_RD_pt[13];
            auto it = AK3_RD_MC5_RD.find(AK3_RD_pt[15]);
            if (it != AK3_RD_MC5_RD.end())
            {
                it->second.push_back(MC5_RD_0_7_10_13);
            }
            else
            {
                vector<array<cell_t, 4> > MC5_RD_0_7_10_13_list = {MC5_RD_0_7_10_13};
                AK3_RD_MC5_RD.insert(pair<cell_t, vector<array<cell_t, 4> > >(AK3_RD_pt[15], MC5_RD_0_7_10_13_list));
            }
        }
    }
    for (auto it = AK3_RD_MC5_RD.begin(); it != AK3_RD_MC5_RD.end(); it++)
    {
        cout << hex << setfill('0') << setw(2) << ((*it).first) + '\0' << "    " << dec << setfill(' ') << ((*it).second.size()) << endl;
    }
    cout << "Real DoFB = 2^" << log2((double)counter[(((uint32_t)(AK3_C_pt[ 0])) | ((uint32_t)(AK3_C_pt[ 5]) << 8) | ((uint32_t)(AK3_C_pt[10]) << 16))]) << endl;
    cout << "For other possible constants, degree of freedom are: " << endl;
    for (size_t i = 0ULL; i < (1ULL<<24); i++)
    {
        cout << hex << setfill('0') << setw(4) << i << "    2^" << dec << setfill(' ') << log2((double)counter[i]) << endl;
    }
}

void testRD()
{
    state_t SB4_RD  = _mm_setzero_si128();
    state_t k3_RD   = _mm_setzero_si128();
    state_t AK3_RD  = _mm_setzero_si128();
    state_t SB5_RD  = _mm_setzero_si128();
    state_t SB6_RD  = _mm_setzero_si128();
    state_t k4_RD   = _mm_setzero_si128();
    state_t k5_RD   = _mm_setzero_si128();

    cell_t * k3_RD_pt   = (cell_t * )&k3_RD  ;
    cell_t * AK3_RD_pt  = (cell_t * )&AK3_RD ;

    for (size_t i = 0; i < SB_a; i++)
    {
        AK3_RD_pt[15] = i;
        k3_RD_pt[15] = i;

        SB4_RD = _mm_set1_epi8(0x52);
         k4_RD = enc_one_round(k3_RD, RC[4]);
        SB5_RD = enc_one_round(SB4_RD, k4_RD);

        k5_RD = enc_one_round(k4_RD, RC[5]);
        SB6_RD = enc_one_round(SB5_RD, k5_RD);
        printstate(SB6_RD);
    }
}

/**************************************************************************************
 * Precompute values of Blue cells (neutral bytes for forward) in AK3,
 * the values all have constant impact on Red cells in MC2, 
 * the values all have constant impact on Red/Gray cells in SB4.
 **************************************************************************************/
vector<array<cell_t, 8> > AK3_BL_list;
array<cell_t, 8> AK3_BL_8cell;
array<array<cell_t, 2>, (1UL << CELL_InBits)> MC3_BL_8_5_list;
array<array<cell_t, 2>, (1UL << CELL_InBits)> MC3_BL_9_6_list;
array<array<cell_t, 2>, (1UL << CELL_InBits)> MC3_BL_4_11_list;
array<array<cell_t, 2>, (1UL << CELL_InBits)> MC3_BL_10_7_list;
array<vector<array<cell_t, 4> >, (1UL << (2UL * CELL_InBits))> AK3_BL_5_10_to_MC3_BL_8_5_9_6_list;
array<vector<array<cell_t, 4> >, (1UL << (2UL * CELL_InBits))> AK3_BL_5_10_to_MC3_BL_4_11_10_7_list;
array<cell_t, 4> mc_in;
array<cell_t, 4> MC3_BL_4cell;

void precomputeBL()
{
    for (uint32_t AK3_BL_1_2_i = 0; AK3_BL_1_2_i < (1UL << (DoFF - 16UL)); AK3_BL_1_2_i++)
    {
        state_t AK3_BL;
        state_t MC3_BL;
        state_t SB3_BL;
        cell_t * AK3_BL_pt = (cell_t *)&AK3_BL;
        cell_t * MC3_BL_pt = (cell_t *)&MC3_BL;
        cell_t * SB3_BL_pt = (cell_t *)&SB3_BL;

        cell_t SB3_BL_8;
        cell_t SB3_BL_9;
        cell_t MC3_BL_8;
        cell_t MC3_BL_5;
        cell_t SB3_BL_13;
        cell_t SB3_BL_14;
        cell_t MC3_BL_9;
        cell_t MC3_BL_6;
        cell_t SB3_BL_4;
        cell_t SB3_BL_7;
        cell_t MC3_BL_4;
        cell_t MC3_BL_11;
        cell_t SB3_BL_2;
        cell_t SB3_BL_3;
        cell_t MC3_BL_10;
        cell_t MC3_BL_7;

        state_t MC3_BL_tmp;
        state_t AK3_BL_tmp;
        cell_t * MC3_BL_tmp_pt = (cell_t *)&MC3_BL_tmp;
        cell_t * AK3_BL_tmp_pt = (cell_t *)&AK3_BL_tmp;

        cell2_t AK3_BL_5_10;

        state_t AK3_BL_determined;
        cell_t * AK3_BL_determined_pt = (cell_t *)&AK3_BL_determined;

        state_t AK3_BL_verify;
        state_t SB3_BL_verify;
        state_t MC2_BL_verify;
        state_t AK3_BL_org;
        state_t MC2_BL_org;

        AK3_BL = AK3_C;
        AK3_BL_pt[1] = AK3_BL_1_2_i & 0xFF;
        AK3_BL_pt[2] = (AK3_BL_1_2_i >> CELL_InBits) & 0xFF;
        MC3_BL = inv_mix_column(AK3_BL);
        SB3_BL = inv_SB_SR(MC3_BL);

        for (uint32_t SB3_BL_8_i = 0; SB3_BL_8_i < (1UL << CELL_InBits); SB3_BL_8_i++)
        { 
            SB3_BL_8 = (cell_t)SB3_BL_8_i;
            mc_in = {SB3_BL_8, SB3_BL_pt[10], SB3_C_pt[11], MC2_C_pt[9]};
            SB3_BL_9 = mix_column_I_4671_O_5(mc_in);
            MC3_BL_8 = S[SB3_BL_8];
            MC3_BL_5 = S[SB3_BL_9];
            MC3_BL_8_5_list[SB3_BL_8_i] = {MC3_BL_8, MC3_BL_5};
        }
        for (uint32_t SB3_BL_13_i = 0; SB3_BL_13_i < (1UL << CELL_InBits); SB3_BL_13_i++)
        { 
            SB3_BL_13 = (cell_t)SB3_BL_13_i;
            mc_in = {SB3_C_pt[12], SB3_BL_13, SB3_BL_pt[15], MC2_C_pt[12]};
            SB3_BL_14 = mix_column_I_4570_O_6(mc_in);
            MC3_BL_9  = S[SB3_BL_13];
            MC3_BL_6  = S[SB3_BL_14];
            MC3_BL_9_6_list[SB3_BL_13_i] = {MC3_BL_9, MC3_BL_6};
        }
        for (uint32_t SB3_BL_4_i = 0; SB3_BL_4_i < (1UL << CELL_InBits); SB3_BL_4_i++)
        { 
            SB3_BL_4 = (cell_t)SB3_BL_4_i;
            mc_in = {SB3_BL_4, SB3_BL_pt[5], SB3_C_pt[6], MC2_C_pt[6]};
            SB3_BL_7  = mix_column_I_4562_O_7(mc_in);
            MC3_BL_4  = S[SB3_BL_4];
            MC3_BL_11 = S[SB3_BL_7];
            MC3_BL_4_11_list[SB3_BL_4_i] = {MC3_BL_4, MC3_BL_11};
        }
        for (uint32_t SB3_BL_2_i = 0; SB3_BL_2_i < (1UL << CELL_InBits); SB3_BL_2_i++)
        { 
            SB3_BL_2 = (cell_t)SB3_BL_2_i;
            mc_in = {SB3_BL_pt[0], SB3_C_pt[1], SB3_BL_2, MC2_C_pt[3]};
            SB3_BL_3  = mix_column_I_4563_O_7(mc_in);
            MC3_BL_10 = S[SB3_BL_2];
            MC3_BL_7  = S[SB3_BL_3];
            MC3_BL_10_7_list[SB3_BL_2_i] = {MC3_BL_10, MC3_BL_7};
        }
        for (auto & MC3_BL_8_5_e: MC3_BL_8_5_list)
        {
            for (auto & MC3_BL_9_6_e: MC3_BL_9_6_list)
            {
                MC3_BL_tmp = _mm_setzero_si128();
                AK3_BL_tmp = _mm_setzero_si128();
                MC3_BL_tmp_pt[8] = MC3_BL_8_5_e[0];
                MC3_BL_tmp_pt[5] = MC3_BL_8_5_e[1];
                MC3_BL_tmp_pt[9] = MC3_BL_9_6_e[0];
                MC3_BL_tmp_pt[6] = MC3_BL_9_6_e[1];
                AK3_BL_tmp = mix_column(MC3_BL_tmp);
                AK3_BL_tmp = xor_states(AK3_BL_tmp, AK3_C);
                AK3_BL_5_10 = ((cell2_t)AK3_BL_tmp_pt[5]) | (((cell2_t)AK3_BL_tmp_pt[10]) << CELL_InBits);
                MC3_BL_4cell = {
                    MC3_BL_tmp_pt[8],
                    MC3_BL_tmp_pt[5],
                    MC3_BL_tmp_pt[9],
                    MC3_BL_tmp_pt[6]
                };
                AK3_BL_5_10_to_MC3_BL_8_5_9_6_list[AK3_BL_5_10].push_back(MC3_BL_4cell);
            }
        }
        for (auto & MC3_BL_4_11_e: MC3_BL_4_11_list)
        {
            for (auto & MC3_BL_10_7_e: MC3_BL_10_7_list)
            {
                MC3_BL_tmp = _mm_setzero_si128();
                AK3_BL_tmp = _mm_setzero_si128();
                MC3_BL_tmp_pt[ 4] = MC3_BL_4_11_e[0];
                MC3_BL_tmp_pt[11] = MC3_BL_4_11_e[1];
                MC3_BL_tmp_pt[10] = MC3_BL_10_7_e[0];
                MC3_BL_tmp_pt[ 7] = MC3_BL_10_7_e[1];
                AK3_BL_tmp = mix_column(MC3_BL_tmp);
                AK3_BL_5_10 = ((cell2_t)AK3_BL_tmp_pt[5]) | (((cell2_t)AK3_BL_tmp_pt[10]) << CELL_InBits);
                MC3_BL_4cell = {
                    MC3_BL_tmp_pt[ 4],
                    MC3_BL_tmp_pt[11],
                    MC3_BL_tmp_pt[10],
                    MC3_BL_tmp_pt[ 7]
                };
                AK3_BL_5_10_to_MC3_BL_4_11_10_7_list[AK3_BL_5_10].push_back(MC3_BL_4cell);
            }
        }
        for (uint32_t AK3_BL_5_10_i = 0; AK3_BL_5_10_i < (1UL << (2UL * CELL_InBits)); AK3_BL_5_10_i++)
        {
            cell2_t AK3_BL_5_10 = (cell2_t)AK3_BL_5_10_i;
            if ((AK3_BL_5_10_to_MC3_BL_8_5_9_6_list[AK3_BL_5_10].size() != 0) && 
                (AK3_BL_5_10_to_MC3_BL_4_11_10_7_list[AK3_BL_5_10].size() != 0))
            {
                for (auto & MC3_BL_8_5_9_6_e : AK3_BL_5_10_to_MC3_BL_8_5_9_6_list[AK3_BL_5_10])
                {
                    for (auto & MC3_BL_4_11_10_7_e : AK3_BL_5_10_to_MC3_BL_4_11_10_7_list[AK3_BL_5_10])
                    {
                        MC3_BL_pt[ 8] = MC3_BL_8_5_9_6_e[0];
                        MC3_BL_pt[ 5] = MC3_BL_8_5_9_6_e[1];
                        MC3_BL_pt[ 9] = MC3_BL_8_5_9_6_e[2];
                        MC3_BL_pt[ 6] = MC3_BL_8_5_9_6_e[3];
                        MC3_BL_pt[ 4] = MC3_BL_4_11_10_7_e[0];
                        MC3_BL_pt[11] = MC3_BL_4_11_10_7_e[1];
                        MC3_BL_pt[10] = MC3_BL_4_11_10_7_e[2];
                        MC3_BL_pt[ 7] = MC3_BL_4_11_10_7_e[3];

                        AK3_BL_determined = mix_column(MC3_BL);
                        AK3_BL_8cell[0] = AK3_BL_determined_pt[1];
                        AK3_BL_8cell[1] = AK3_BL_determined_pt[2];
                        AK3_BL_8cell[2] = AK3_BL_determined_pt[4];
                        AK3_BL_8cell[3] = AK3_BL_determined_pt[6];
                        AK3_BL_8cell[4] = AK3_BL_determined_pt[7];
                        AK3_BL_8cell[5] = AK3_BL_determined_pt[8];
                        AK3_BL_8cell[6] = AK3_BL_determined_pt[9];
                        AK3_BL_8cell[7] = AK3_BL_determined_pt[11];
                        AK3_BL_list.push_back(AK3_BL_8cell);
                        #if DEBUG
                        AK3_BL_verify = mix_column(MC3_BL);
                        AK3_BL_verify = and_states(AK3_BL_verify, AK3_C_MASK);
                        AK3_BL_org    = and_states(AK3_BL, AK3_C_MASK);
                        if (!EQU(AK3_BL_verify, AK3_BL_org)) cout << "MC3_BL to AK3_BL wrong." << endl;
                        SB3_BL_verify = inv_SB_SR(MC3_BL);
                        SB3_BL_verify = and_states(SB3_BL_verify, SB3_BL_MASK);
                        SB3_BL_verify = xor_states(SB3_BL_verify, SB3_C);
                        MC2_BL_verify = inv_mix_column(SB3_BL_verify);
                        MC2_BL_verify = and_states(MC2_BL_verify, MC2_C_MASK);
                        MC2_BL_org    = and_states(MC2_C, MC2_C_MASK);
                        if (!EQU(MC2_BL_verify, MC2_BL_org)) cout << "MC3_BL to MC2_C wrong." << endl;
                        #endif
                    }
                }
            }
            AK3_BL_5_10_to_MC3_BL_8_5_9_6_list[AK3_BL_5_10].clear();
            AK3_BL_5_10_to_MC3_BL_4_11_10_7_list[AK3_BL_5_10].clear();
        }
    }
    cout << "2^" << log2((double)(AK3_BL_list.size())) << endl;

    #if DEBUG
    cout << hex << setfill('0');
    array<unsigned char, 4> in;
    unsigned char O1;
    unsigned char O2;
    for (size_t i = 0ULL; i < (1ULL<<32); i++)
    {
        in[0] = (i >> (0 * 8)) & 0xFFU;
        in[1] = (i >> (1 * 8)) & 0xFFU;
        in[2] = (i >> (2 * 8)) & 0xFFU;
        in[3] = (i >> (3 * 8)) & 0xFFU;
        O1 = mix_column_I_4671_O_5(in);
        O2 = mix_column_I_4671_O_5_test(in);
        if (O1 != O2) cout << setw(2) << O1 + '\0' << " != "  << setw(2) << O2 + '\0' << endl;
        O1 = mix_column_I_4570_O_6(in);
        O2 = mix_column_I_4570_O_6_test(in);
        if (O1 != O2) cout << setw(2) << O1 + '\0' << " != "  << setw(2) << O2 + '\0' << endl;
        O1 = mix_column_I_4563_O_7(in);
        O2 = mix_column_I_4563_O_7_test(in);
        if (O1 != O2) cout << setw(2) << O1 + '\0' << " != "  << setw(2) << O2 + '\0' << endl;
        O1 = mix_column_I_4562_O_7(in);
        O2 = mix_column_I_4562_O_7_test(in);
        if (O1 != O2) cout << setw(2) << O1 + '\0' << " != "  << setw(2) << O2 + '\0' << endl;
        O1 = mix_column_I_0126_O_3(in);
        O2 = mix_column_I_0126_O_3_test(in);
        if (O1 != O2) cout << setw(2) << O1 + '\0' << " != "  << setw(2) << O2 + '\0' << endl;
        O1 = mix_column_I_1235_O_0(in);
        O2 = mix_column_I_1235_O_0_test(in);
        if (O1 != O2) cout << setw(2) << O1 + '\0' << " != "  << setw(2) << O2 + '\0' << endl;
    }
    #endif
}

vector<uint32_t>* Forward_L[omp_nb_threads];

/**************************************************************************************
 * Attack
 **************************************************************************************/
bool attack(state_t Target)
{
    clock_t timecnt;

    uint64_t complexity[omp_nb_threads] = {0ULL};
    uint64_t complexity_F[omp_nb_threads] = {0ULL};
    uint64_t complexity_B[omp_nb_threads] = {0ULL};
    uint64_t complexity_M[omp_nb_threads] = {0ULL};

    uint64_t complexity_Sum = 0ULL;
    uint64_t complexity_F_Sum = 0ULL;
    uint64_t complexity_B_Sum = 0ULL;
    uint64_t complexity_M_Sum = 0ULL;

    timecnt = clock();
    precomputeBL();
    timecnt = clock() - timecnt;
    cout << "Precompute neutral bytes for forward takes " << (((float)timecnt)/(CLOCKS_PER_SEC * 60))  << " (mins)" << endl;
    cout << endl;

    bool findFlag = false;
    timecnt = clock();
    for (int ti = 0; ti < omp_nb_threads; ti++)
    {
        Forward_L[ti] = new vector<uint32_t>[1ULL<<32ULL];
    }
    #pragma omp parallel for num_threads(omp_nb_threads)
    for (uint64_t k3_C_i = 0;  k3_C_i < CONST_SAMPLE_N; k3_C_i++)  // sampling for the CONST_SAMPLE_N Gray bytes in k3_C
    {
        // private variables for each thread
        int tid = omp_get_thread_num();
        clock_t timecnt_thisthread;
        fout[tid].open(fn[tid], ios::app);

        // Fixed constant in key state
        state_t k3_C;
        state_t k4_C;
        state_t KMC4_C;
        cell_t * k3_C_pt = (cell_t *)&k3_C; // pointer to Constants in #k3_C

        // variables for forward (Blue + Gray) independently computing
        state_t AK3_BL;
        state_t SB4_BL;
        state_t SB5_BL;
        state_t MC5_BL;
        state_t SB6_BL;
        state_t MC6_BL;
        state_t AK6_BL;
        cell_t * AK3_BL_pt = (cell_t *)&AK3_BL;
        cell_t * MC5_BL_pt = (cell_t *)&MC5_BL;
        cell_t * AK6_BL_pt = (cell_t *)&AK6_BL;
        cell4_t BL_match;

        // variables for backward (Red + Gray) independently computing
        state_t km_RD = _mm_setzero_si128();
        state_t k0_RD = _mm_setzero_si128();
        state_t k1_RD = _mm_setzero_si128();
        state_t k2_RD = _mm_setzero_si128();
        state_t k3_RD = _mm_setzero_si128();
        state_t k4_RD = _mm_setzero_si128();
        state_t k5_RD = _mm_setzero_si128();
        state_t k6_RD = _mm_setzero_si128();
        cell_t * k3_RD_pt  = (cell_t * )&k3_RD;
        state_t AK3_RD;
        state_t SB3_RD;
        state_t AK2_RD;
        state_t MC2_RD;
        state_t SB2_RD;
        state_t AK1_RD;
        state_t MC1_RD;
        state_t SB1_RD;
        state_t AK0_RD;
        state_t AK0_RD_Guessed;
        state_t MC0_RD;
        state_t SB0_RD;
        state_t AK6_RD;
        cell_t * AK0_RD_pt = (cell_t *)&AK0_RD;
        cell_t * AK3_RD_pt = (cell_t *)&AK3_RD;
        cell_t * AK6_RD_pt = (cell_t *)&AK6_RD;
        cell4_t RD_match;
        size_t BL_RD_matched_n;
        size_t AK3_BL_idx;

        // variables for combined computing after match at the matching point
        uint32_t matchedcell = 0;
        state_t AK3_BL_MT;
        state_t MC3_BL_MT;
        state_t SB3_BL_MT;
        state_t SB3_RD_MT;
        state_t SB3;
        state_t SB2;
        state_t SB1;
        state_t AK0;
        state_t AK0_Guessed;
        state_t SB0;
        state_t AK6_Backward;
        state_t SB4;
        state_t SB5;
        state_t SB6;
        state_t AK6_Forward;
        cell_t * AK3_BL_MT_pt = (cell_t *)&AK3_BL_MT;

        // variables for final check after finding the target
        state_t SStart_check;
        state_t KStart_check;
        state_t k3_check;
        state_t k4_check;
        state_t k5_check;
        state_t k6_check;
        state_t k2_check;
        state_t k1_check;
        state_t k0_check;
        state_t km_check;
        state_t State_BWD_check;
        state_t State_FWD_check;

        k3_C = k3_fixed_C; // Constants in k3
        //Fix k3[0, 5, 10, 15] to be iS[0] ^ AK3[0, 5, 10, 15]
        k3_C_pt[ 1] = (k3_C_i >> (0 * 8)) & 0xFF;
        k3_C_pt[ 2] = (k3_C_i >> (1 * 8)) & 0xFF;
        k3_C_pt[ 4] = (k3_C_i >> (2 * 8)) & 0xFF;
        k3_C_pt[ 6] = (k3_C_i >> (3 * 8)) & 0xFF;
        k3_C_pt[ 7] = (k3_C_i >> (4 * 8)) & 0xFF;
        k3_C_pt[ 8] = (k3_C_i >> (5 * 8)) & 0xFF;
        k3_C_pt[ 9] = (k3_C_i >> (6 * 8)) & 0xFF;
        k3_C_pt[11] = (k3_C_i >> (7 * 8)) & 0xFF;
        k4_C   = enc_one_round(k3_C, RC[4]);
        KMC4_C = SB_SR(k4_C);

        // Forward compute -- Blue cells
        #if PRINT_COUNT
        fout[tid] << "tid " << tid << ": Forward compute -- Blue cells" << endl;
        #endif
        for (uint64_t Forward_L_i = 0ULL; Forward_L_i < (1ULL << 32ULL); Forward_L_i++)
        {
            Forward_L[tid][Forward_L_i].clear();
        }
        uint64_t AK3_BL_n = AK3_BL_list.size();
        for (uint64_t AK3_BL_i = 0ULL; AK3_BL_i < AK3_BL_n; AK3_BL_i++)
        {
            AK3_BL = AK3_C;
            AK3_BL_pt[ 1] = AK3_BL_list[AK3_BL_i][0];
            AK3_BL_pt[ 2] = AK3_BL_list[AK3_BL_i][1];
            AK3_BL_pt[ 4] = AK3_BL_list[AK3_BL_i][2];
            AK3_BL_pt[ 6] = AK3_BL_list[AK3_BL_i][3];
            AK3_BL_pt[ 7] = AK3_BL_list[AK3_BL_i][4];
            AK3_BL_pt[ 8] = AK3_BL_list[AK3_BL_i][5];
            AK3_BL_pt[ 9] = AK3_BL_list[AK3_BL_i][6];
            AK3_BL_pt[11] = AK3_BL_list[AK3_BL_i][7];

            SB4_BL = add_round_key(AK3_BL, k3_C);
            SB5_BL = enc_one_round(SB4_BL, k4_C);
            MC5_BL = SB_SR(SB5_BL);
            MC5_BL = xor_states(MC5_BL, KMC4_C);
            #if DEBUG
            if ((MC5_BL_pt[0] != 0) || (MC5_BL_pt[7] != 0) || (MC5_BL_pt[10] != 0) || (MC5_BL_pt[13] != 0))
            {
                fout[tid] << "Forward compute to MC5_BL wrong." << endl;
                fout[tid].close();
                exit(1);
            }
            #endif
            SB6_BL = mix_column(MC5_BL);
            SB6_BL = xor_states(SB6_BL, RC[5]);
            MC6_BL = SB_SR(SB6_BL);
            AK6_BL = mix_column(MC6_BL);
            BL_match = 
                (((cell4_t)(AK6_BL_pt[ 0]) & 0xFFUL) << (0UL * CELL_InBits)) |
                (((cell4_t)(AK6_BL_pt[ 5]) & 0xFFUL) << (1UL * CELL_InBits)) |
                (((cell4_t)(AK6_BL_pt[10]) & 0xFFUL) << (2UL * CELL_InBits)) |
                (((cell4_t)(AK6_BL_pt[15]) & 0xFFUL) << (3UL * CELL_InBits));
            Forward_L[tid][BL_match].push_back((uint32_t)AK3_BL_i);
            complexity_F[tid]++;
        }
        #if PRINT_COUNT
        fout[tid] << "tid " << tid << ": Forward compute -- Done" << endl << endl;
        #endif

        // Backward compute -- Red cells
        #if PRINT_COUNT
        fout[tid] << "tid " << tid << ": Backward compute -- Red cells" << endl;
        #endif
        for (uint64_t k3_RD_15_init = 0ULL; k3_RD_15_init < (1ULL<<DoFB); k3_RD_15_init++)
        {
            AK3_RD        = AK3_C;
            k3_RD         = k3_C;
            k3_RD_pt[15]  ^= k3_RD_15_init;
            AK3_RD_pt[15] ^= k3_RD_15_init;

            k4_RD = enc_one_round(k3_RD, RC[4]);
            k5_RD = enc_one_round(k4_RD, RC[5]);
            k6_RD = enc_one_round(k5_RD, RC[6]);
            k2_RD = dec_one_round(k3_RD, RC[3]);
            k1_RD = dec_one_round(k2_RD, RC[2]);
            k0_RD = dec_one_round(k1_RD, RC[1]);
            #if (MP_Mode==0)
            km_RD = dec_one_round(k0_RD, RC[0]);
            #endif


            SB3_RD = dec_one_round(AK3_RD, zero_128);
            SB3_RD = and_states(SB3_RD, SB3_RD_MASK);
            AK2_RD = add_round_key(SB3_RD, k2_RD);
            MC2_RD = inv_mix_column(AK2_RD);
            MC2_RD = xor_states(MC2_RD, MC2_C);
            SB2_RD = inv_SB_SR(MC2_RD);
            SB1_RD = dec_one_round(SB2_RD, k1_RD);
            AK0_RD = add_round_key(SB1_RD, k0_RD);
            AK0_RD = and_states(AK0_RD, AK0_RD_MASK);
            for (uint64_t AK0_RD_guessed_i = 0ULL; AK0_RD_guessed_i < (1ULL<<GUESS_RD_InBits); AK0_RD_guessed_i++)
            {
                AK0_RD_pt[0] = (AK0_RD_guessed_i >> (0UL * CELL_InBits)) & 0xFF;
                AK0_RD_pt[2] = (AK0_RD_guessed_i >> (1UL * CELL_InBits)) & 0xFF;
                AK0_RD_pt[3] = (AK0_RD_guessed_i >> (2UL * CELL_InBits)) & 0xFF;
                AK0_RD_Guessed = and_states(AK0_RD, AK0_GUESS_MASK);
                MC0_RD = inv_mix_column(AK0_RD);
                SB0_RD = inv_SB_SR(MC0_RD);
                #if (MP_Mode==0)
                AK6_RD = add_round_key(SB0_RD, km_RD);
                #endif
                AK6_RD = add_round_key(AK6_RD, Target);
                AK6_RD = add_round_key(AK6_RD, k6_RD);
                
                RD_match = 
                (((cell4_t)(AK6_RD_pt[ 0]) & 0xFFUL) << (0UL * CELL_InBits)) |
                (((cell4_t)(AK6_RD_pt[ 5]) & 0xFFUL) << (1UL * CELL_InBits)) |
                (((cell4_t)(AK6_RD_pt[10]) & 0xFFUL) << (2UL * CELL_InBits)) |
                (((cell4_t)(AK6_RD_pt[15]) & 0xFFUL) << (3UL * CELL_InBits));
                complexity_B[tid]++;

                matchedcell = 0;
                BL_RD_matched_n = Forward_L[tid][RD_match].size();
                if (BL_RD_matched_n != 0)
                {
                    //fout[tid] << "tid " << tid << ": Find match at matching point." << endl;
                    matchedcell += BL_RD_matched_n;
                    for (size_t BL_RD_matched_i = 0; BL_RD_matched_i < BL_RD_matched_n; BL_RD_matched_i++)
                    {
                        complexity_M[tid]++;
                        AK3_BL_MT = AK3_C;
                        AK3_BL_idx = Forward_L[tid][RD_match][BL_RD_matched_i];
                        AK3_BL_MT_pt[ 1] = AK3_BL_list[AK3_BL_idx][0];
                        AK3_BL_MT_pt[ 2] = AK3_BL_list[AK3_BL_idx][1];
                        AK3_BL_MT_pt[ 4] = AK3_BL_list[AK3_BL_idx][2];
                        AK3_BL_MT_pt[ 6] = AK3_BL_list[AK3_BL_idx][3];
                        AK3_BL_MT_pt[ 7] = AK3_BL_list[AK3_BL_idx][4];
                        AK3_BL_MT_pt[ 8] = AK3_BL_list[AK3_BL_idx][5];
                        AK3_BL_MT_pt[ 9] = AK3_BL_list[AK3_BL_idx][6];
                        AK3_BL_MT_pt[11] = AK3_BL_list[AK3_BL_idx][7];
                        AK3_BL_MT_pt[15] ^= k3_RD_15_init;
                        MC3_BL_MT = inv_mix_column(AK3_BL_MT);
                        SB3_BL_MT = inv_SB_SR(MC3_BL_MT);
                        SB3_BL_MT = and_states(SB3_BL_MT, SB3_BL_MASK);
                        SB3_RD_MT = and_states(SB3_RD, SB3_RD_MASK);
                        SB3 = xor_states(SB3_BL_MT, SB3_RD_MT);
                        SB2 = dec_one_round(SB3, k2_RD);
                        SB1 = dec_one_round(SB2, k1_RD);
                        AK0 = add_round_key(SB1, k0_RD);
                        AK0_Guessed = and_states(AK0, AK0_GUESS_MASK);
                        if (!EQU(AK0_Guessed, AK0_RD_Guessed))
                        {
                            continue;
                        }
                        //fout[tid] << "tid " << tid << ": Find match at guessing point." << endl;
                        SB0 = dec_one_round(SB1, k0_RD);
                        #if (MP_Mode==0)
                        AK6_Backward = add_round_key(SB0, km_RD);
                        #endif
                        AK6_Backward = add_round_key(AK6_Backward, k6_RD);
                        AK6_Backward = add_round_key(AK6_Backward, Target);
                        AK6_Backward = and_states(AK6_Backward, PARTIAL_MATCH_TARGET_MASK);

                        SB4 = enc_one_round(SB3, k3_RD);
                        SB5 = enc_one_round(SB4, k4_RD);
                        SB6 = enc_one_round(SB5, k5_RD);
                        AK6_Forward = enc_one_round(SB6, zero_128);
                        AK6_Forward = and_states(AK6_Forward, PARTIAL_MATCH_TARGET_MASK);

                        if (EQU(AK6_Backward, AK6_Forward))
                        {
                            #pragma omp critical
                            {
                                cout << "tid " << tid << ": Find match at partial target." << endl;
                                findFlag = true;
                                for (int ti = 0; ti < omp_nb_threads; ti++)
                                {
                                    complexity_F_Sum += complexity_F[ti];
                                    complexity_B_Sum += complexity_B[ti];
                                    complexity_M_Sum += complexity_M[ti];
                                }
                                cout << "Find match on " << PARTIAL_InBits << " bits." << endl;
                                cout << "Complexity (Total Forward ) is 2^" << log2((double)complexity_F_Sum) << endl;
                                cout << "Complexity (Total Backward) is 2^" << log2((double)complexity_B_Sum) << endl;
                                cout << "Complexity (After Matched ) is 2^" << log2((double)complexity_M_Sum) << endl;
                                cout << "Complexity (Total) is 2^" << log2(((double)complexity_F_Sum + (double)complexity_B_Sum)/2.0 + (double)complexity_M_Sum) << endl;
                                timecnt_thisthread = clock() - timecnt;
                                cout << "Find match on " << PARTIAL_InBits << " bits takes " << (((float)timecnt_thisthread)/(CLOCKS_PER_SEC * 60))  << " (mins)" << endl;
                                cout << endl;

                                SStart_check = SB3;
                                KStart_check = k3_RD;
                                k3_check = KStart_check;
                                k4_check = enc_one_round(k3_check, RC[4]);
                                k5_check = enc_one_round(k4_check, RC[5]);
                                k6_check = enc_one_round(k5_check, RC[6]);
                                k2_check = dec_one_round(k3_check, RC[3]);
                                k1_check = dec_one_round(k2_check, RC[2]);
                                k0_check = dec_one_round(k1_check, RC[1]);
                                #if (MP_Mode==0)
                                km_check = dec_one_round(k0_check, RC[0]);
                                #endif

                                cout << "==== Backward computation ====" << endl;
                                State_BWD_check = SStart_check;
                                cout << "#SB3 = " << endl;
                                printstate(State_BWD_check);
                                cout << "k2 = " << endl;
                                printstate(k2_check);
                                State_BWD_check = dec_one_round(State_BWD_check, k2_check);
                                cout << "#SB2 = " << endl;
                                printstate(State_BWD_check);
                                cout << "k1 = " << endl;
                                printstate(k1_check);
                                State_BWD_check = dec_one_round(State_BWD_check, k1_check);
                                cout << "#SB1 = " << endl;
                                printstate(State_BWD_check);
                                cout << "k0 = " << endl;
                                printstate(k0_check);
                                State_BWD_check = dec_one_round(State_BWD_check, k0_check);
                                cout << "#SB0 = " << endl;
                                printstate(State_BWD_check);
                                #if (MP_Mode==0)
                                cout << "km = " << endl;
                                printstate(km_check);
                                State_BWD_check = add_round_key(State_BWD_check, km_check);
                                cout << "k6 = " << endl;
                                printstate(k6_check);
                                State_BWD_check = add_round_key(State_BWD_check, k6_check);
                                cout << "#SB0 ⊕ km ⊕ k6 = " << endl;
                                printstate(State_BWD_check);
                                #else
                                cout << "k6 = " << endl;
                                printstate(k6_check);
                                State_BWD_check = add_round_key(State_BWD_check, k6_check);
                                cout << "#SB0 ⊕ k6 = " << endl;
                                printstate(State_BWD_check);
                                #endif


                                cout << "==== Forward computation ====" << endl;
                                State_FWD_check = SStart_check;
                                cout << "#SB3 = " << endl;
                                printstate(State_FWD_check);
                                cout << "k3 = " << endl;
                                printstate(k3_check);
                                State_FWD_check = enc_one_round(State_FWD_check, k3_check);
                                cout << "#SB4 = " << endl;
                                printstate(State_FWD_check);
                                cout << "k4 = " << endl;
                                printstate(k4_check);
                                State_FWD_check = enc_one_round(State_FWD_check, k4_check);
                                cout << "#SB5 = " << endl;
                                printstate(State_FWD_check);
                                cout << "k5 = " << endl;
                                printstate(k5_check);
                                State_FWD_check = enc_one_round(State_FWD_check, k5_check);
                                cout << "#SB6 = " << endl;
                                printstate(State_FWD_check);
                                State_FWD_check = enc_one_round(State_FWD_check, zero_128);
                                cout << "#AK6 = " << endl;
                                printstate(State_FWD_check);

                                State_FWD_check = xor_states(State_FWD_check, State_BWD_check);
                                #if (MP_Mode==0)
                                cout << "#AK6 ⊕ #SB0 ⊕ km ⊕ k6 = " << endl;
                                #else
                                cout << "#AK6 ⊕ #SB0 ⊕ k6 = " << endl;
                                #endif
                                printstate(State_FWD_check);

                                cout << "Target = " << endl;
                                printstate(Target);

                                fout[tid].close();
                                exit(0);
                            }
                        }
                    }
                }
                #pragma omp critical
                if (findFlag)
                {
                    fout[tid].close();
                    exit(0);
                }
            }
            #pragma omp critical
            if (findFlag)
            {
                fout[tid].close();
                exit(0);
            }
        }
        #pragma omp critical
        if (findFlag)
        {
            fout[tid].close();
            exit(0);
        }
        for (uint64_t Forward_L_i = 0ULL; Forward_L_i < (1ULL << 32ULL); Forward_L_i++)
        {
            Forward_L[tid][Forward_L_i].clear();
        }
        fout[tid].close();
    }
    if (!findFlag)
    {
        for (int ti = 0; ti < omp_nb_threads; ti++)
        {
            complexity_F_Sum += complexity_F[ti];
            complexity_B_Sum += complexity_B[ti];
            complexity_M_Sum += complexity_M[ti];
        }
        cout << "Did not find match on " << PARTIAL_InBits << " bits." << endl;
        cout << "Complexity (Total Forward ) is 2^" << log2((double)complexity_F_Sum) << endl;
        cout << "Complexity (Total Backward) is 2^" << log2((double)complexity_B_Sum) << endl;
        cout << "Complexity (After Matched ) is 2^" << log2((double)complexity_M_Sum) << endl;
        cout << "Complexity (Total) is 2^" << log2(((double)complexity_F_Sum + (double)complexity_B_Sum)/2.0 + (double)complexity_M_Sum) << endl;
        timecnt = clock() - timecnt;
        cout << "Find match on " << PARTIAL_InBits << " bits takes " << (((float)timecnt)/(CLOCKS_PER_SEC * 60))  << " (mins)" << endl;
        cout << endl;
    }
    for (int ti = 0; ti < omp_nb_threads; ti++)
    {
        delete[] Forward_L[ti];
    }
    return findFlag;
}

int main()
{
    for (int ti = 0; ti < omp_nb_threads; ti++)
    {
        fn[ti] = to_string(PARTIAL_InBits) + "_tid_" + to_string(ti) + ".txt";
    }
    //testRD1();
    //testRD();

    // Here, the value of the target is set to be 0s.
    // One can set it to be arbitrary values using block_rndfill:
    state_t T = _mm_setzero_si128();
    // block_rndfill((uint8_t *)&T, sizeof(state_t));

    attack(T);

    return 0;
}