#pragma once
#ifndef SHORT_CHARACTERISTICS_CALCULATIONS_H
#define SHORT_CHARACTERISTICS_CALCULATIONS_H

#include "struct_short_characteristics_headers.h"
#include "struct_short_characteristics_global_structure.h"


size_t ReadFileVtk(const size_t class_file_vtk, const std::string name_file_vtk, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid,
	vtkDataArray*& density, vtkDataArray*& absorp_coef, vtkDataArray*& rad_en_loose_rate, const bool is_print = false);

size_t ReadSphereDirectionDecartToSpherical(const std::string name_file_sphere_direction, vector<Vector3>& directions_all, vector<Type>& squares, Type& square_surface);

int InitGlobalValue(Vector3& start_point_plane_coord, Matrix3& transform_matrix, Matrix3& inverse_transform_matrix,
	Matrix3& straight_face, Matrix3& inclined_face);

int FindNeighborsPairFaceAndBoundaries(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<int>& all_pairs_face);
int FindNeighborsPairFace(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<int>& all_pairs_face);
int GetNumberNeighborFace(const int a, const int b, const int c, vtkCell* neighbor_cell);

int InitNodesValue(const std::vector<int>& all_pairs_face, std::vector<cell>& nodes_value);
int ResetNodesValue(std::vector<cell>& nodes_value);
int ReadNormalFile(std::string& name_file_normals, std::vector<Normals>& normals);

Type NormIllum(const std::vector<Type>& Illum, const std::vector<Type>& Illum2);
int ReadGraph(const std::string name_file_graph, std::vector<int>& sorted_id_cell);
int ReadGraphBin(const std::string name_file_graph, std::vector<int>& sorted_id_cell);
int SetVertexMatrix(const size_t number_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, Eigen::Matrix4d& vertex_tetra);
int FindInAndOutFaces(const Vector3& direction, const int number_cell, const std::vector<Normals>& normals, int* face_state);
size_t CenterOfTetra(const int number_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Vector3& point_in_tetra);
size_t FindAllCenterOfTetra(const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, std::vector<Vector3>& point_in_tetra);

int FromGlobalToLocalTetra(const Eigen::Matrix4d& vertex_tetra, const Eigen::Vector3d& global_coord, Eigen::Vector3d& local_coord);
int FromLocalToGlobalTetra(const Eigen::Matrix4d& vertex_tetra, const Eigen::Vector3d& local_coord, Eigen::Vector3d& global_coord);

int FromTetraToPlane(const Eigen::Matrix3d& transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& tetra_coord, Eigen::Vector3d& plane_coord);
int FromPlaneToTetra(const Eigen::Matrix3d& inverse_transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& plane_coord,
	Eigen::Vector3d& tetra_coord);

size_t IntersectionWithPlane(vtkCell* face, const Vector3& start_point, const Vector3& direction, Vector3& result);
bool InTriangle(const int num_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cell_face, int number_face, const Eigen::Vector3d& XX);

Type BoundaryFunction(const int id_cell, const Vector3& x, const Vector3& direction, const std::vector<Type>& illum_old,
	const vector<Vector3>& directions, const vector<Type>& squares);

Vector3 GetInterpolationCoef(const Eigen::Matrix3d& interpolation_nodes, const Eigen::Vector3d& function_value);

Type GetS(const int num_cell, const Vector3& direction, const std::vector<Type>& illum_old,
	const vector<Vector3>& directions, const vector<Type>& squares);
int Min(const int a, const int b);

size_t MakeEnergy(const vector<Type>& Illum, const vector<Type>& squares, const Type scuare_surface, vector<Type>& energy);
Type IntegarteDirection(const int num_cell, const vector<Type>& Illum, const vector<Type>& squares, const Type scuare_surface);

size_t WriteFileSolution(const std::string name_file_out, const std::vector<Type>& vector_illum, const std::vector<Type>& vector_energy,
	vtkSmartPointer<vtkUnstructuredGrid>& u_grid);
#endif