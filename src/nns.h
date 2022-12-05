/*
 * nns.h
 *
 *  Created on: 6 oct. 2022
 *      Author: Kévin Carrier
 */

#ifndef NNS_H_
#define NNS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// red_cost_model : model for NNS (from Martin Albretch code)
// These are not asymptotic expressions but compress the data in [AC:AGPS20]_ with the fix and improvement from [MATZOV22]_ applied which covers up to beta = 1024
double NN_AGPS[21][2] = {
		{0.4215069316732438, 20.166968300536567}, // all_pairs-classical
		{0.3171724396445733, 25.2982895173379}, // all_pairs-dw
		{0.31552858350028, 22.478746811528104}, // all_pairs-g
		{0.3222895263943547, 36.11746438609664}, // all_pairs-ge19
		{0.41862512941897706, 9.899382685790897}, // all_pairs-naive_classical
		{0.31401512571180035, 7.694659414353819}, // all_pairs-naive_quantum
		{0.31553282513562797, 20.87859415484879}, // all_pairs-t_count
		{0.29613500308205365, 20.387885985467914}, // list_decoding-classical
		{0.2663676536352464, 25.299541499216627}, // list_decoding-dw
		{0.26600114174341505, 23.440974518186337}, // list_decoding-g
		{0.26799889622667994, 30.839871638418543}, // list_decoding-ge19
		{0.29371310617068064, 15.930690682515422}, // list_decoding-naive_classical
		{0.2632557273632713, 15.685687713591548}, // list_decoding-naive_quantum
		{0.2660264010780807, 22.432158856991474}, // list_decoding-t_count
		{0.3558614423344473, 23.08252781663665}, // random_buckets-classical
		{0.30704199602260734, 25.58196897625173}, // random_buckets-dw
		{0.30610964725102396, 22.928235564044588}, // random_buckets-g
		{0.31089687605567917, 36.02129974535213}, // random_buckets-ge19
		{0.35448283789554536, 15.28878540793911}, // random_buckets-naive_classical
		{0.3021142178390157, 11.151745066682524}, // random_buckets-naive_quantum
		{0.3061477007403873, 21.418301489775203}, // random_buckets-t_count
};

#define all_pairs_classical 0
#define all_pairs_dw 1
#define all_pairs_g 2
#define all_pairs_ge19 3
#define all_pairs_naive_classical 4
#define all_pairs_naive_quantum 5
#define all_pairs_t_count 6
#define list_decoding_classical 7
#define list_decoding_dw 8
#define list_decoding_g 9
#define list_decoding_ge19 10
#define list_decoding_naive_classical 11
#define list_decoding_naive_quantum 12
#define list_decoding_t_count 13
#define random_buckets_classical 14
#define random_buckets_dw 15
#define random_buckets_g 16
#define random_buckets_ge19 17
#define random_buckets_naive_classical 18
#define random_buckets_naive_quantum 19
#define random_buckets_t_count 20


/*
 *  @brief Nearest-Neigbhor Search cost
 *  @param beta dimension
 *  @param red_cost_model model
 *  @return The time complexity of the  Nearest-Neigbhor Search.
 */
double Tnns_f(const double beta, const int red_cost_model){
	return pow(2.0, NN_AGPS[red_cost_model][0] * beta + NN_AGPS[red_cost_model][1]);
}

#endif /* NNS_H_ */
