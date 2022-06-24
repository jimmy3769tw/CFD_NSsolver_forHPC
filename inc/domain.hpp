#pragma once
#include <vector>
#include <tuple>
#include <iostream>
class calDomain{

	public:
	// ----------------------------
	int word_size;

	int i_begin, j_begin, k_begin, i_endof, j_endof, k_endof;

	std::vector<int> 
	i_begin_table, 	j_begin_table, 	k_begin_table, 
	i_endof_table, 	j_endof_table, 	k_endof_table, 
	i_length_table,	j_length_table, k_length_table;
	// ----------------------------

	void initLocal(int i, int j, int k);

	void show_begin();

	void initTable( std::vector<int> gridASize, std::vector<int> dims);

	std::pair<std::vector<int>, std::vector<int> >
	get_begEnd(int slid){
		auto beg = i_begin_table;
		auto end = i_endof_table;
		// -------------------
		for(auto &x: beg){ x -=2;}
		for(auto &x: end){ x -=2;}
		for(auto &x: beg){ x *=slid;}
		for(auto &x: end){ x *=slid;}
		// -------------------
		return std::make_pair(beg, end);
	}

	private:
	int gC = 2;
	void divideDomain( int dims, int n, vector<int> &begin ,vector<int> &endof ,vector<int> &length );

};


void calDomain::initLocal(int i, int j, int k){
	i_begin = i_begin_table.at(i);
	j_begin = j_begin_table.at(j);
	k_begin = k_begin_table.at(k);

	i_endof =i_endof_table.at(i);
	j_endof =j_endof_table.at(j);
	k_endof =k_endof_table.at(k);
}



void calDomain::initTable(	
	std::vector<int> gridASize, 
	std::vector<int> dims
){
	word_size = dims.at(0) * dims.at(1) * dims.at(2);
	divideDomain(dims.at(0), gridASize.at(0)-2*gC, i_begin_table, i_endof_table, i_length_table);
	divideDomain(dims.at(1), gridASize.at(1)-2*gC, j_begin_table, j_endof_table, j_length_table);
	divideDomain(dims.at(2), gridASize.at(2)-2*gC, k_begin_table, k_endof_table, k_length_table);
}


void calDomain::show_begin(){
	cout << "show begin";

	for (auto x :i_begin_table){
		std::cout << x << ", ";
	}

	for (auto x :j_begin_table){
		std::cout << x << ", ";
	}

	for (auto x :k_begin_table){
		std::cout << x << ", ";
	}
}


void calDomain::divideDomain(
	int dims, int n,
	std::vector<int> &begin ,
	std::vector<int> &endof ,
	std::vector<int> &length
)
{
	if (dims == 0){
		throw std::invalid_argument( "divideDomain !!!!\n");
	}

	begin.resize(dims);

	endof.resize(dims);
	
	length.resize(dims);

	for (int i = 0; i < dims; ++i)
	{
		length[i] = (n - n % dims) / (dims);
	}

	for (int i = 0; i < n % dims; ++i)
	{
		length[i] += 1;
	}

	int gC = 2;

	begin[0] = gC;

	for (int i = 0; i < dims - 1; i++)
	{
		begin[i + 1] = begin[i] + length[i];
	}

	for (int i = 0; i < dims; i++)
	{
		endof[i] = begin[i] + length[i];
	}

}
