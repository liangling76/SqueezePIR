#include <iostream>
#include <unistd.h>
#include <cmath>
#include <chrono>
#include "seal/seal.h"


#define POLY_DEGREE 8192
#define NUM_SLOT 4096

#define NUM_OBJ 40960
#define NUM_COL 4096     
#define RANK 32          
#define IDX 2069

#define NUM_BLOCK NUM_OBJ / NUM_SLOT

#define CHECK_CORRECT true


typedef struct SqueezePirDB{
    std::vector<std::vector<std::vector<double>>> A; // (NUM_BLOCK, R, SLOT)
    std::vector<std::vector<std::vector<double>>> B; // (NUM_BLOCK, R, SLOT)
}SqueezePirDB;


// generate the input query 
std::vector<seal::Ciphertext> query_gen(seal::CKKSEncoder &encoder, seal::Encryptor &encryptor, double scale){
    std::vector<seal::Ciphertext> query(NUM_BLOCK);

    std::vector<double> tmp_msg(NUM_SLOT);
    seal::Plaintext tmp_pt;

    // encrypt query 
    for(int i = 0; i < NUM_BLOCK; i++)
    {
        tmp_msg = std::vector<double>(NUM_SLOT, 0ULL);
        for(int j = 0; j < NUM_SLOT; j++)
        {
            if(i * NUM_SLOT + j == IDX) tmp_msg[j] = 1;
        }
        encoder.encode(tmp_msg, scale, tmp_pt);
        encryptor.encrypt_symmetric(tmp_pt, query[i]);
    }

    
    return query;
}


// generate the database for SqueezePIR with matrix A and B
SqueezePirDB squeezepir_db_gen(){
    auto rng = std::default_random_engine {};

    SqueezePirDB db;
    db.A = std::vector<std::vector<std::vector<double>>>(NUM_BLOCK);
    db.B = std::vector<std::vector<std::vector<double>>>(NUM_BLOCK);

    for(int i = 0; i < NUM_BLOCK; i++)
    {
        db.A[i] = std::vector<std::vector<double>>(RANK, std::vector<double>(NUM_SLOT, 0ULL));
        db.B[i] = std::vector<std::vector<double>>(RANK, std::vector<double>(NUM_SLOT, 0ULL));
    
        // db.A
        for(int j = 0; j < NUM_SLOT; j++)
        {
            double gamma_sqr_remain = 1.0;

            for(int k = 0; k < RANK; k++)
            {
                double gamma = ((double) rand() / double(RAND_MAX)) * sqrt(gamma_sqr_remain);
                if((double) rand() / double(RAND_MAX) > 0.5) gamma *= -1;                    
                db.A[i][k][j] = gamma;
                gamma_sqr_remain -= gamma * gamma;
            }
        }

        // db.B
        for(int j = 0; j < NUM_COL; j++)
        {
            double gamma_sqr_remain = 1.0;

            for(int k = 0; k < RANK; k++)
            {
                double gamma = ((double) rand() / double(RAND_MAX)) * sqrt(gamma_sqr_remain);
                if((double) rand() / double(RAND_MAX) > 0.5) gamma *= -1;
                db.B[i][k][j] = gamma;
                gamma_sqr_remain -= gamma * gamma;
            }
        }
    
    }

    return db;
}



int main(){

    // std::cout << seal::CoeffModulus::MaxBitCount(16384) << "\n";

    /* settings for ckks*/
    seal::EncryptionParameters parms(seal::scheme_type::ckks);
    parms.set_poly_modulus_degree(POLY_DEGREE);
    parms.set_coeff_modulus(seal::CoeffModulus::Create(POLY_DEGREE, { 60, 49, 49, 60 }));
    double scale = pow(2.0, 49);

    seal::SEALContext context(parms);

    seal::KeyGenerator keygen(context);
    seal::PublicKey public_key;
    seal::RelinKeys relin_keys;
    seal::GaloisKeys gal_keys;

    auto secret_key = keygen.secret_key();
    keygen.create_relin_keys(relin_keys);
    keygen.create_galois_keys(gal_keys);

    seal::Encryptor encryptor(context, secret_key);
    seal::Evaluator evaluator(context);
    seal::Decryptor decryptor(context, secret_key);

    seal::CKKSEncoder encoder(context);


    /* time initial */
    std::chrono::high_resolution_clock::time_point time_start_db, time_end_db, time_start_squeezepir, time_end_squeezepir;  
    std::chrono::high_resolution_clock::time_point time_expand1, time_qa_dba, time_expand2, time_qb_dbb;  



    /* build the query */
    std::vector<seal::Ciphertext> query = query_gen(encoder, encryptor, scale);


    /* build the database */
    SqueezePirDB db = squeezepir_db_gen();

    time_start_db = std::chrono::high_resolution_clock::now();
    
    std::vector<std::vector<seal::Ciphertext>> db_A(NUM_BLOCK, std::vector<seal::Ciphertext>(RANK));
    std::vector<std::vector<seal::Ciphertext>> db_B(NUM_BLOCK, std::vector<seal::Ciphertext>(RANK));

    seal::Plaintext tmp_pt;

    for(int i = 0; i < NUM_BLOCK; i++)
    {
        for(int k = 0; k < RANK; k++)
        {
            encoder.encode(db.A[i][k], scale, tmp_pt);
            encryptor.encrypt_symmetric(tmp_pt, db_A[i][k]);

            encoder.encode(db.B[i][k], scale, tmp_pt);
            encryptor.encrypt_symmetric(tmp_pt, db_B[i][k]);
        }
    }

    time_end_db = std::chrono::high_resolution_clock::now();
    auto tmp_time = (std::chrono::duration_cast<std::chrono::microseconds>(time_end_db - time_start_db)).count();
    std::cout << "data processing time: " << tmp_time * 1e-6<< "\n";


    /* data retrieving */
    time_start_squeezepir = std::chrono::high_resolution_clock::now();

    seal::Ciphertext tmp_ct_A, tmp_ct_B, tmp_ct;

    int expand_times = (int)log2(NUM_SLOT);

    
    // comptue query_exp = expand(query)
    std::vector<seal::Ciphertext> query_exp(NUM_BLOCK);
    
    for(int i = 0; i < NUM_BLOCK; i++)
    {
        evaluator.rotate_vector(query[i], 1, gal_keys, tmp_ct);
        evaluator.add(query[i], tmp_ct, query_exp[i]);
        for(int j = 1; j < expand_times; j++)
        {
            int rotate_step = (int)pow(2.0, j);
            evaluator.rotate_vector(query_exp[i], rotate_step, gal_keys, tmp_ct);
            evaluator.add_inplace(query_exp[i], tmp_ct);
        }
    }

    time_expand1 = std::chrono::high_resolution_clock::now();
    tmp_time = (std::chrono::duration_cast<std::chrono::microseconds>(time_expand1 - time_start_squeezepir)).count();
    std::cout << "expand query time: " << tmp_time  * 1e-6 << "\n";


    // compute result_A = sum(query * db_A)
    std::vector<seal::Ciphertext> result_A(RANK);
    std::vector<bool> tmp_inplace_A(RANK, false);

    for(int i = 0; i < NUM_BLOCK; i++)
    {
        for(int k = 0; k < RANK; k++)
        {
            if(!tmp_inplace_A[k])
            {
                evaluator.multiply(query[i], db_A[i][k], result_A[k]);
                // evaluator.relinearize_inplace(result_A[k], relin_keys);
                tmp_inplace_A[k] = true;
            }
            else
            {
                evaluator.multiply(query[i], db_A[i][k], tmp_ct_A);
                // evaluator.relinearize_inplace(tmp_ct_A, relin_keys);
                evaluator.add_inplace(result_A[k], tmp_ct_A);
            }
        }
    }

    for(int k = 0; k < RANK; k++)
    {
        evaluator.relinearize_inplace(result_A[k], relin_keys);
        evaluator.rescale_to_next_inplace(result_A[k]);
    }

    time_qa_dba = std::chrono::high_resolution_clock::now();
    tmp_time = (std::chrono::duration_cast<std::chrono::microseconds>(time_qa_dba - time_expand1)).count();
    // std::cout << "q_a db_a time: " << tmp_time << "\n";


    // compute result_B = sum(query_exp * db_B)
    std::vector<seal::Ciphertext> result_B(RANK);
    std::vector<bool> tmp_inplace_B(RANK, false);

    for(int i = 0; i < NUM_BLOCK; i++)
    {
        for(int k = 0; k < RANK; k++)
        {
            if(!tmp_inplace_B[k])
            {
                evaluator.multiply(query_exp[i], db_B[i][k], result_B[k]);
                // evaluator.relinearize_inplace(result_B[k], relin_keys);
                tmp_inplace_B[k] = true;
            }
            else
            {
                evaluator.multiply(query_exp[i], db_B[i][k], tmp_ct);
                // evaluator.relinearize_inplace(tmp_ct, relin_keys);
                evaluator.add_inplace(result_B[k], tmp_ct);
            }
        }
    }

    for(int k = 0; k < RANK; k++)
    {
        evaluator.relinearize_inplace(result_B[k], relin_keys);
        evaluator.rescale_to_next_inplace(result_B[k]);
    }

    time_qb_dbb = std::chrono::high_resolution_clock::now();
    tmp_time = (std::chrono::duration_cast<std::chrono::microseconds>(time_qb_dbb - time_expand1)).count();
    // std::cout << "q_b db_b time: " << tmp_time << "\n";
    std::cout << "qa qb time: " << tmp_time  * 1e-6 << "\n";


    // comptue result_A = expand(result_A)
    for(int k = 0; k < RANK; k++)
    {
        for(int i = 0; i < expand_times; i++){
            int rotate_step = (int)pow(2.0, i);
            evaluator.rotate_vector(result_A[k], rotate_step, gal_keys, tmp_ct);
            evaluator.add_inplace(result_A[k], tmp_ct);
        }
    }

    time_expand2 = std::chrono::high_resolution_clock::now();
    tmp_time = (std::chrono::duration_cast<std::chrono::microseconds>(time_expand2 - time_qb_dbb)).count();
    std::cout << "expand result A: " << tmp_time  * 1e-6 << "\n";


    // compute result = result_A * result_B
    seal::Ciphertext result;

    for(int k = 0; k < RANK; k++)
    {
        if(k == 0)
        {
            evaluator.multiply(result_A[k], result_B[k], result);
            evaluator.relinearize_inplace(result, relin_keys);
        }
        else
        {
            evaluator.multiply(result_A[k], result_B[k], tmp_ct);
            evaluator.relinearize_inplace(tmp_ct, relin_keys);
            evaluator.add_inplace(result, tmp_ct);
        }
    }
    // evaluator.relinearize_inplace(result, relin_keys);
    evaluator.rescale_to_next_inplace(result);

    time_end_squeezepir = std::chrono::high_resolution_clock::now();
    tmp_time = (std::chrono::duration_cast<std::chrono::microseconds>(time_end_squeezepir - time_expand2)).count();
    std::cout << "result time: " << tmp_time  * 1e-6 << "\n\n";


    auto squeezepir_time = (std::chrono::duration_cast<std::chrono::microseconds>(time_end_squeezepir - time_start_squeezepir)).count();
    std::cout << "\nSqueezePIR time: " << squeezepir_time * 1e-6 << "\n";


    /* check correctness */
    if(CHECK_CORRECT)
    {
        std::cout << "decryption...\n";
        // decrypt
        std::vector<double> result_msg;
        decryptor.decrypt(result, tmp_pt);
        encoder.decode(tmp_pt, result_msg);

        int db_A_i = IDX / NUM_SLOT;
        int db_A_j = IDX % NUM_SLOT;
        int db_B_i = IDX / NUM_SLOT;

        int err_4 = 0, err_5 = 0, err_6 = 0, err_7 = 0, err_8 = 0, err_9 = 0;

        for(int j = 0; j < NUM_COL; j++)
        {
            double org = 0;
            double msg = result_msg[j];

            for(int k = 0; k < RANK; k++)
            {
                org += db.A[db_A_i][k][db_A_j] * db.B[db_B_i][k][j];
            }

        
            if (org - msg > 1e-4 || msg - org > 1e-4) err_4 ++;

            if (org - msg > 1e-5 || msg - org > 1e-5) err_5 ++;

            if (org - msg > 1e-6 || msg - org > 1e-6) err_6 ++;

            if (org - msg > 1e-7 || msg - org > 1e-7) err_7 ++;

            if (org - msg > 1e-8 || msg - org > 1e-8) err_8 ++;

            if (org - msg > 1e-9 || msg - org > 1e-9) err_9 ++;

            // std::cout << j << "\t" << org << "\t" << msg << "\n";
        }

        std::cout << "error larger than 1e-4: " << err_4 << "\n";
        std::cout << "error larger than 1e-5: " << err_5 << "\n";
        std::cout << "error larger than 1e-6: " << err_6 << "\n";
        std::cout << "error larger than 1e-7: " << err_7 << "\n";
        std::cout << "error larger than 1e-8: " << err_8 << "\n";
        std::cout << "error larger than 1e-9: " << err_9 << "\n";
        
    }



    return 0;
}

