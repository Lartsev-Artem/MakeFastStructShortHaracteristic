#pragma once
#ifndef SHORT_CHARACTERISTICS_LOGIC_H
#define SHORT_CHARACTERISTICS_LOGIC_H

#include "struct_short_characteristics_global_structure.h"
#include "struct_short_characteristics_headers.h"
#include "struct_short_characteristics_calculations.h"

size_t IntersectionWithPlaneDisk(const Vector3& X0, const Vector3& n, Vector3& res);
int GetNodes(const int num_cur_cell, const std::vector<Face>& grid, const ShortId num_cur_out_face,
	const Eigen::Matrix4d& vertex_tetra, const int* face_state, const Vector3& direction, const std::vector<Normals>& normals,
	std::vector<cell>& nodes_value,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_res_bound,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_s,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_x,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_x0_local,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_in_id);
int GetNodes(const int num_cur_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell, const int num_cur_out_face,
	const Eigen::Matrix4d& vertex_tetra, const int* face_state, const Vector3& direction,
	std::vector<cell>& nodes_value, const std::vector<Type>& illum_old,
	const vector<Vector3>& directions, const vector<Type>& squares);
int GetNodes(const int num_cur_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell, const ShortId num_cur_out_face,
	const Eigen::Matrix4d& vertex_tetra, const int* face_state, const Vector3& direction,
	std::vector<cell>& nodes_value,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_res_bound,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_s,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_x,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_x0_local,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_in_id);

int CalculateNodeValue(const int num_cur_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell,
	const int num_cur_out_face, const int* face_state,
	const Vector3& direction, std::vector<cell>& nodes_value,
	const int num_node, const Eigen::Matrix4d& vertex_tetra, Vector3& x,
	const std::vector<Type>& illum_old, const vector<Vector3>& directions, const vector<Type>& squares);

int CalculateNodeValue(const int num_cur_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell,
	const int num_cur_out_face, const int* face_state, const Vector3& direction,
	std::vector<cell>& nodes_value, const int num_node, const Eigen::Matrix4d& vertex_tetra, Vector3& x,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_res_bound,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_s,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_x0_local,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_in_id);

Type CalculateIllumeOnInnerFace(const int num_cell, const int num_in_face, const Eigen::Matrix4d& vertex_tetra,
	const Vector3& x, const Vector3& x0, const std::vector<cell>& nodes_value);
Type CurGetIllum(const int cur_id, const Vector3 x, const Type s, const Type I_node_prev, const Vector3& cur_direction,
	const std::vector<Type>& illum_old,
	const vector<Vector3>& directions, const vector<Type>& squares);

Type CalculateIllumeOnInnerFace(const int num_cell, const int num_in_face, const Eigen::Matrix4d& vertex_tetra,
	const Vector3& x, const Vector3& x0, const std::vector<cell>& nodes_value,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_res_bound,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_x0_local);

Type GetValueInCenterCell(const int num_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell, const Vector3 center,
	const Vector3 direction,
	const Eigen::Matrix4d& vertex_tetra, const std::vector<cell>& nodes_value,
	const std::vector<Type>& illum_old, const vector<Vector3>& directions, const vector<Type>& squares);

#endif
