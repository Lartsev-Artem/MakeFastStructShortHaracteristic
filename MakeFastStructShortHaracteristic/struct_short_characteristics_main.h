#pragma once
#ifndef SHORT_CHARACTERISTICS_MAIN_H
#define SHORT_CHARACTERISTICS_MAIN_H

#include "struct_short_characteristics_headers.h"
#include "struct_short_characteristics_global_structure.h"
#include "struct_short_characteristics_calculations.h"
#include "struct_short_characteristics_logic_function.h"

template<typename Type>
size_t ReadStartSettings(Type name_file_settings, std::string& name_file_vtk,
	std::string& name_file_sphere_direction, std::string& out_file_grid_vtk, std::string& main_dir) {

	std::ifstream ifile;
	ifile.open(name_file_settings);
	if (!ifile.is_open()) {
		std::cerr << " Error : file settings build graph is not open !\n";
		return 1;
	}

	std::string str; // переменная для перевода строки при чтении из файла

	getline(ifile, name_file_vtk);
	getline(ifile, name_file_sphere_direction);
	getline(ifile, out_file_grid_vtk);
	getline(ifile, main_dir);

	ifile.close();
	return 0;
}

template <typename T>
int ReadBinGeneralFile(const std::string& name_file, const std::vector<T>& data) {

	FILE* f;
	f = fopen(name_file.c_str(), "rb");
	if (!f) RETURN("ReadBinGeneralFile not open\n")
		int n;
	fread_unlocked(&n, sizeof(int), 1, f);
	data.resize(n);
	fread_unlocked(data.data(), sizeof(T), n, f);
	fclose(f);

	return 0;
}
#endif