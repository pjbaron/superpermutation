// SuperPermutation.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <stack>


#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))


/// best for 7 = 5906
const int N = 6;

/// how many cuts per chromosone during crossover
/// NOTE: must be less than N
const int splitPoints = 3;

/// how many individuals can there be over and beyond the essential permutation set
const int populationGrowth = 1024;

/// decrease fitness of older individuals
const double eugenics = 0.001;

/// limit the chromosomal length to the calculated lower bounds * this multiplier
const double lengthCapMultiplier = 2.0;

/// percentage chance of picking a permutation as a parent instead of using another evolved chromosome
/// (the permutations set contains every combination of symbols, so this is a form of mutation)
const double usePermutationParent = 0.05;

/// percentage chance of mutating a single location on a chromosome
/// TODO: use a 'wave' which periodically increases this then decreases it again to stimulate exploration
///const double mutateOneLocationChanceStart = 0.05;


class List;
class Individual;

// global variables
int generation = 0;
List** permutations;
Individual** population;
int populationSize;
double mutateOneLocationChance = 0;



static unsigned long rx = 123456789, ry = 362436069, rz = 521288629;
unsigned long ul_rnd()
{
	unsigned long t;
	rx ^= rx << 16;
	rx ^= rx >> 5;
	rx ^= rx << 1;

	t = rx;
	rx = ry;
	ry = rz;
	rz = t ^ rx ^ ry;

	return rz;
}

double rnd()
{
	double v = ul_rnd();
	return v / ((double)ULONG_MAX + 1.0);
}




class List
{
public:
	int length;
	int* data;
public:
	List() { length = 0; data = nullptr; }
	List(int _length) { length = _length; data = new int[length]; }
	List(int _length, int* _data) { length = _length; data = new int[length]; clone(_data); }
	List(List* _list) { length = _list->length; data = new int[length]; clone(_list->data); }
	~List() { if (data != nullptr) delete[] data; data = nullptr; length = 0; }

	void clone(int* _data)
	{
		for (int i = 0; i < length; i++)
			data[i] = _data[i];
	}

	int* clone()
	{
		int* ret = new int[length];
		for (int i = 0; i < length; i++)
			ret[i] = data[i];
		return ret;
	}

	int* slice(int _first, int _length)
	{
		int last = _first + _length;
		_ASSERT(last <= length);
		int* ret = new int[_length];
		for (int j = 0, i = _first; i < last; j++, i++)
		{
			ret[j] = data[i];
		}
		return ret;
	}

	void splice(int _first, int _length)
	{
		int newLength = length - _length;
		int last = _first + _length;
		_ASSERT(last <= length);

		int* ret = new int[newLength];

		int i, j;
		for (i = 0; i < _first; i++)
			ret[i] = data[i];
		for (j = last; j < length; j++, i++)
			ret[i] = data[j];

		delete[] data;
		data = ret;
		length = newLength;
	}

	List* concat(List* _list)
	{
		int i, l = length + _list->length;
		List* ret = new List(l);

		for (i = 0; i < length; i++)
			ret->data[i] = data[i];

		for (int j = 0; i < l; i++, j++)
			ret->data[i] = _list->data[j];

		return ret;
	}

	void swapData(int _i, int _j)
	{
		int t = data[_i];
		data[_i] = data[_j];
		data[_j] = t;
	}
};


class Individual
{
public:
	double fitness;
	int permutations;
	int birthday;
	List* data;
public:
	Individual() { fitness = 0; permutations = 0; birthday = 0; data = nullptr; }
	Individual(double _fitness, List* _list, int _birthday) { fitness = _fitness; permutations = 0; birthday = _birthday; data = new List(_list);  }
	~Individual() { delete data; data = nullptr; fitness = 0; permutations = 0; }
};



// create one of every permutation of symbols in _n into global 'permutations'
int heapPermutation(List* _list, int _size, int _n, int _p)
{
	if (_size == 1)
	{
		int* s = _list->slice(0, _n);
		permutations[_p] = new List(_n, s);
		delete [] s;
		return _p + 1;
	}

	for (int i = 0; i < _size; i++)
	{
		_p = heapPermutation(_list, _size - 1, _n, _p);
		if (_size & 1)
			_list->swapData(0, _size - 1);
		else
			_list->swapData(i, _size - 1);
	}
	return _p;
}

//void integerValue(_list)
//{
//	int t = 0;
//	for (int i = 0, l = _list.length; i < l; i++)
//		t += _list[i] * Math.pow(10, i);
//	return t;
//}

int factorial(int _n)
{
	int t = 1;
	for (int i = 1; i <= _n; i++)
		t *= i;
	return t;
}

int lowestBound(int _n)
{
	int t = factorial(_n);
	t += factorial(_n - 1);
	t += factorial(_n - 2);
	t += factorial(_n - 3);
	t += _n - 3;
	return t;
}

double* rouletteWheel(Individual** _population, int popSize, double& total)
{
	double* _roulette = new double[popSize];
	total = 0;
	for (int i = 0; i < popSize; i++)
	{
		total += _population[i]->fitness;
		_roulette[i] = total;
	}
	return _roulette;
}

Individual* pickParent(Individual** _population, int _popSize, double* _roulette, double _rouletteTotal, int _numPermutations, bool _usePermutation)
{
	// percentage chance that we'll pick a raw permutation for this parent
	if (_usePermutation)
	{
		int i = (int)floor(rnd() * (double)_numPermutations);
		return _population[i];
	}

	// otherwise use the roulette wheel to pick one from the whole set (including raw permutations)
	int i = 0;
	double j = rnd() * _rouletteTotal;

	// https://en.wikipedia.org/wiki/Binary_search_algorithm
	int l = 0, r = _popSize - 1;
	while (r != l) 
	{
		int m = (l + r + 1) / 2;	// https://stackoverflow.com/questions/2745074/fast-ceiling-of-an-integer-division-in-c-c
		if (_roulette[m] > j)
			r = m - 1;
		else
			l = m;
	}
	return _population[l];
	//while (j > _roulette[i])
	//	i++;
	//return _population[i];
}

List* crossOver(List* _parent1, List* _parent2)
{
	int* child = new int[_parent1->length + _parent2->length];
	int* splits1 = new int[splitPoints];
	int* splits2 = new int[splitPoints];
	int j1 = 0, j2 = 0, c = 0;
	for (int i = 0; i < splitPoints; i++)
	{
		// find split points in both parents
		int p1, p2;
		int j1 = 0, j2 = 0;
		if (i > 0) j1 = splits1[i - 1];
		if (i > 0) j2 = splits2[i - 1];
		do {
			p1 = (int)floor(rnd() * (_parent1->length - (splitPoints - i - 1))) + 1;
		} while (p1 <= j1);
		do {
			p2 = (int)floor(rnd() * (_parent2->length - (splitPoints - i - 1))) + 1;
		} while (p2 <= j2);

		// store them so we can grab the previous split point on the next loop
		splits1[i] = p1;
		splits2[i] = p2;

		if (rnd() < 0.5)
		{
			// take from parent1 to p1
			for (; j1 < p1; j1++)
				child[c++] = _parent1->data[j1];
		}
		else
		{
			// take from parent2 to p2
			for (; j2 < p2; j2++)
				child[c++] = _parent2->data[j2];
		}
	}

	// clone the first c elements of child into a new List object to return
	List* ret = new List(c, child);

	delete[] splits2;
	delete[] splits1;
	delete[] child;

	return ret;
}

Individual* createChild(Individual* _parent1, Individual* _parent2)
{
	List* xover = crossOver(_parent1->data, _parent2->data);
	Individual* child = new Individual(0, xover, generation);
	delete xover;
	return child;
}

Individual* createChildByAdding(Individual* _parent1, Individual* _parent2)
{
	List* c = _parent1->data->concat(_parent2->data);
	Individual* child = new Individual(0, c, generation);
	delete c;
	return child;
}

// look for _value in _list searching from _start to _last (not end-inclusive)
int findValue(int _value, int* _list, int _start, int _last)
{
	for (int i = _start; i < _last; i++)
		if (_list[i] == _value)
			return i;
	return -1;
}

// squash them by removing overlapping sequences (1,2,3,2,3,4 => 1,2,3,4) up to length N:
// for each value in the list with i
//  for each matching value in the list up to N positions away with j
//   advance i & j by one
//   if we have advanced i == original j then we have an overlapping sequence
//   else compare values, if they are not equal, break inner loop
void squash(List* _list)
{
	for (int i = 0, l = _list->length - N; i < l; i++)
	{
		int ii = i;

		bool overlap = false;
		int j = i;
		while ((j = findValue(_list->data[i], _list->data, j + 1, i + N)) != -1)
		{
			int jj = j;

			do {
				ii++;
				jj++;
				if (_list->data[ii] != _list->data[jj])
				{
					// no overlap here
					break;
				}
			} while (ii < j);
			if (ii >= j)
			{
				overlap = true;
				// break the inner loop, we're removing i..j so search no more inside that span
				break;
			}
		}

		if (overlap)
		{
			// we found an overlap from i to j
			_list->splice(i, j - i);
			l = _list->length - N;
		}
	}
}

void mutate(List* _list)
{
	if (rnd() < mutateOneLocationChance)
	{
		if (rnd() < 0.5)
		{
			// mutate location i
			int i = (int)floor(rnd() * _list->length);
			_list->data[i] = (int)floor(rnd() * N + 1);
		}
		else
		{
			// swap location i and j
			int i = (int)floor(rnd() * _list->length);
			int j;
			do {
				j = (int)floor(rnd() * _list->length);
			} while (j == i);
			int t = _list->data[i];
			_list->data[i] = _list->data[j];
			_list->data[j] = t;
		}
	}
}

// count how many entire consecutive _permutations exist in _list
int countPermutations(List** _permutations, List* _list, const int _numPermutations)
{
	int* pointers = new int[_numPermutations] {};
	int c = 0;

	// for every value in this _list
	for (int i = 0, l = _list->length; i < l; i++)
	{
		int v = _list->data[i];

		// for every permutation
		for (int j = 0; j < _numPermutations; j++)
		{
			// if we haven't matched this permutation yet and there's still room for it to fit before the end of _list
			if (pointers[j] < N && i < l - N + 1)
			{
				// if the value matches the next permutation value
				if (v == _permutations[j]->data[pointers[j]])
				{
					// advance the pointer to the next permutation value
					if (++pointers[j] >= N)
					{
						c++;		// complete match, count it
					}
				}
				else
				{
					// reset to the start
					pointers[j] = 0;
					// if the value matches the start permutation value
					if (v == _permutations[j]->data[pointers[j]])
					{
						// advance the pointer to the next permutation value
						pointers[j]++;
					}
				}
			}
		}

		// exit if we've found all of the permutations
		if (c >= _numPermutations)
			break;
	}

	delete[] pointers;
	return c;
}

/*
// count how many entire consecutive _permutations exist in _list
int countPermutations(List** _permutations, List* _list, int _numPermutations)
{
	int c = 0;
	// for every permutation
	for (int i = 0; i < _numPermutations; i++)
	{
		List* permutation = _permutations[i];
		// for every value in this _list
		for (int j = 0, m = _list->length - permutation->length; j < m; j++)
		{
			// for every value in the permutation
			int k, n;
			for (k = 0, n = permutation->length; k < n; k++)
			{
				if (_list->data[j + k] != permutation->data[k])
				{
					// they don't match, gtfo
					break;
				}
			}

			// if we found this permutation
			if (k >= n)
			{
				// count it once, then break to look for the next permutation
				c++;
				break;
			}
		}
	}
	return c;
}
*/

/*
int findWorstIndividual(Individual** _population, int _popSize, int _numPermutations)
{
	int worstIndividual = -1;
	double worstFitness = DBL_MAX;
	for (int i = _numPermutations; i < _popSize; i++)
	{
		// find the one with the worst fitness
		if (_population[i]->fitness < worstFitness)
		{
			worstFitness = _population[i]->fitness;
			worstIndividual = i;
		}

		// if it's a tie, find the shortest one
		if (worstIndividual != -1)
		{
			if (_population[i]->fitness == worstFitness && _population[i]->data->length < _population[worstIndividual]->data->length)
			{
				worstFitness = _population[i]->fitness;
				worstIndividual = i;
			}
		}
	}
	return worstIndividual;
}

int reducePopulation(Individual** _population, int _popSize, int _numPermutations, int _reduceBy)
{
	for (int i = 0; i < _reduceBy; i++)
	{
		int j = findWorstIndividual(_population, _popSize, _numPermutations);
		if (j == -1) break;

		// remove element j from the _population
		delete _population[j];
		for (int k = j; k < _popSize - 1; k++)
			_population[k] = _population[k + 1];
		_popSize--;
	}
	return _popSize;
}
*/

int cullSortedPopulation(Individual** _population, int _popSize, int _numPermutations, int _reduceBy)
{
	int j = _popSize - _reduceBy;
	while(j < _popSize)
	{
		delete _population[j];
		_population[j] = nullptr;
		j++;
	}
	return _popSize - _reduceBy;
}

double populationDiversity(int _numPermutations)
{
	// population diversity measure from two randomly chosen 'pin' locations dropping through the evolved population
	int d = 0;
	int pin1 = (int)floor(rnd() * N);
	int pin2;
	do {
		pin2 = (int)floor(rnd() * N);
	} while (pin1 == pin2);

	for (int i = _numPermutations; i < populationSize - 1; i++)
	{
		if (population[i]->data->data[pin1] != population[i + 1]->data->data[pin1])
			d++;

		if (population[i]->data->length > pin2)
		{
			if (population[i + 1]->data->length > pin2)
			{
				// both valid, compare them
				if (population[i]->data->data[pin2] != population[i + 1]->data->data[pin2])
					d++;
			}
			else
			{
				// first valid, second invalid
				d++;
			}
		}
		else
		{
			// first invalid, second valid
			if (population[i + 1]->data->length > pin2)
				d++;
		}
	}

	return (double)d;
}





// qsort comparitor
int compareFitnesses(const void * a, const void * b)
{
	Individual* ia = *(Individual **)a;
	Individual* ib = *(Individual **)b;
	return (int)(ib->fitness - ia->fitness);
}

int main()
{
/*	// validate the rng
	double minr = DBL_MAX;
	double maxr = DBL_MIN;
	for (int i = 0; i < 1000000; i++)
	{
		double r = rnd();
		if (r < minr) minr = r;
		if (r > maxr) maxr = r;
	}
	std::cout << "random min " << minr << " max " << maxr << std::endl;
*/

	const int numPermutations = factorial(N);
	const int maximumPopulation = (numPermutations + populationGrowth) * 2;

	/// NOTE: good results = n! + (n-1)! + (n-2)! + (n-3!) + n - 3  (Egan et al)
	const double goodResult = lowestBound(N);

	///the maximum incentive for a genome to grow longer
	const double maxIncentive = (goodResult * goodResult);


	permutations = new List*[ numPermutations ];

	// recursive, but maximum stack depth = N
	int symbols[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };
	List* symbolList = new List(20, symbols);
	_ASSERT( N <= 20 );

	heapPermutation(symbolList, N, N, 0);

	std::cout << "permutations of " << N << " = " << numPermutations << std::endl;
	std::cout << "lowest bound = " << goodResult << std::endl;

	// assign one permutation to each member of the initial population
	population = new Individual*[ maximumPopulation ];
	for (int i = 0; i < numPermutations; i++)
	{
		population[i] = new Individual( 0.01, permutations[i], generation );
	}
	populationSize = numPermutations;

	// GA
	// while any result contains all permutations with (length < n! + (n-1)! + (n-2)! + (n-3!) + n - 3) and is short enough
	//  sort population for ascending fitness
	//  cull the population to bring it down to the size limit
	//  for each new child to be created
	//   pick parents using fitness values
	//   combine parents to make a child by slicing them both and mixing the slices
	//   squash the child to remove overlapping sequences
	//   evaluate child fitness by counting the number of permutations
	//  end for
	// end while

	while(populationSize == numPermutations || population[numPermutations]->data->length >= goodResult || population[numPermutations]->permutations < numPermutations)
	{
		qsort(&(population[numPermutations]), populationSize - numPermutations, sizeof(Individual*), compareFitnesses);
		if (populationSize - numPermutations >= populationGrowth)
		{
			int reduceBy = populationSize - numPermutations - populationGrowth + 1;
			populationSize = cullSortedPopulation(population, populationSize, numPermutations, reduceBy);
			int p = population[numPermutations]->permutations;
			int l = population[numPermutations]->data->length;
			int a = generation - population[numPermutations]->birthday;
			if (p == numPermutations) std::cout << "*";
			std::cout << "fittest = " << population[numPermutations]->fitness << " permutations = " << p << " length = " << l << " age = " << a << std::endl;
		}

		if (populationSize > numPermutations + 1)
		{
			double d = populationDiversity( numPermutations );
			mutateOneLocationChance = (100.0 / (double)d);
			mutateOneLocationChance = min(mutateOneLocationChance, 0.5);
			std::cout << "diversity = " << d << " mutation chance = " << mutateOneLocationChance << std::endl;
		}

		if (eugenics > 0)
		{
			// age negatively affects fitness
			for (int i = numPermutations; i < populationSize; i++)
			{
				population[i]->fitness -= eugenics;
				population[i]->fitness = max(population[i]->fitness, 0.01);
			}
		}

		std::cout << "generation = " << generation << " population = " << populationSize << std::endl;

		double total;
		double* roulette = rouletteWheel(population, populationSize, total /* out */);
		int p = populationSize;
		for (int i = 0; i < populationGrowth; i++)
		{
			Individual* parent1 = pickParent(population, populationSize, roulette, total, numPermutations, false);
			Individual* parent2 = pickParent(population, populationSize, roulette, total, numPermutations, rnd() < usePermutationParent);
			Individual* child = createChild(parent1, parent2);
			mutate(child->data);
			squash(child->data);

			// discard results that are too long
			if (child->data->length <= (int)((double)goodResult * lengthCapMultiplier))
			{
				double incentiveToGrow = (double)child->data->length - (double)goodResult;
				double lengthFitness = (maxIncentive - incentiveToGrow * incentiveToGrow) / maxIncentive;
				child->permutations = countPermutations(permutations, child->data, numPermutations);
				child->fitness = (double)child->permutations + lengthFitness;
				population[p++] = child;
			}
			else
			{
				delete child;
			}
		}

		delete [] roulette;
		populationSize = p;
		generation++;
	}

	std::cout << "\nResult!" << std::endl;
	for (int i = 0; i < population[numPermutations]->data->length; i++)
	{
		std::cout << population[numPermutations]->data->data[i] << ",";
	}
	std::cout << std::endl;
	std::cout << "fitness = " << population[numPermutations]->fitness << " length = " << population[numPermutations]->data->length << std::endl;
	std::cout << "done." << std::endl;

	return 0;
}
