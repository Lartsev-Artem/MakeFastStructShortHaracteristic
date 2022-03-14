#include "struct_short_characteristics_main.h"
Vector3 start_point_plane_coord;   // начало координат плоскости
Matrix3 transform_matrix;          // матрица перехода из базового тетраэдра в плоскость
Matrix3 inverse_transform_matrix;  // матрица перехода из плоскости в базовый тетраэдр

Matrix3	straight_face;  // 3 узла интерпол€ции
Matrix3 inclined_face;  // 3 узла интерпол€ции на наклонной плоскости

//std::vector<Type> Illum2;

// скал€рные данные сетки (unstructured_grid)
vtkDataArray* density;
vtkDataArray* absorp_coef;
vtkDataArray* rad_en_loose_rate;

Type square_surface;  // площадь поверхности дискретной 

Vector3 center_local_sphere;  // центр описанной сферы около стандартного тетраэдра

int num_cur_direction; // номер текущего направлени€
Vector3 cur_direction;

int count_negative_interpolation; // число отрицательных значений интерпол€ции интесивности

std::vector<Vector3> X;
std::vector<Vector2> X0;
std::vector<Type> S;

std::vector<Type> res_inner_bound;
std::vector<Vector3> x_try_surface;
std::vector<int> id_try_surface;

std::string BASE_ADRESS;

ShortId n_out;

int pos_x_try;
int posX;
int posX0;
int posOutC;
int posOut;
int posIn;
int posS;
int posRes;



//const Vector3 center_point(0, 0, 0);
//const Type inner_radius =  0.51; // радиус внутренней сферы (с запасом)

const Vector3 center_point(1, 0, 0);
const Type inner_radius = 0.11; // радиус внутренней сферы (с запасом)


int main(int argc, char* argv[])
{
	std::string name_file_settings = "";

	if (argc <= 1)
		name_file_settings = "D:\\Desktop\\FilesCourse\\settings_file_make.txt";
	else
		name_file_settings = argv[1];


	std::string name_file_vtk;
	std::string name_file_sphere_direction;
	std::string name_file_graph;
	std::string out_file_grid_vtk;
	std::string name_file_normals;
	

	if (ReadStartSettings(name_file_settings, name_file_vtk, name_file_sphere_direction, out_file_grid_vtk, BASE_ADRESS, name_file_normals)) {
		std::cout << "Error reading the start settings\n";
		return 1;
	}

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();


	Type _clock = -omp_get_wtime();
	if (ReadFileVtk(0, name_file_vtk, unstructured_grid, density, absorp_coef, rad_en_loose_rate, true)) {
		std::cout << "Error reading the file vtk\n";
		return 1;
	}
	_clock += omp_get_wtime();
	std::cout << "\n Reading time of the vtk_grid file: " << _clock << "\n";

	// make
	std::vector<Face> grid;
	std::vector<Eigen::Matrix4d> vertexs;
	{
		std::string name_file_cells = BASE_ADRESS + "grid.bin";
		std::string name_file_vertex = BASE_ADRESS + "vertex.bin";		
		
		ReadCellFaces(name_file_cells, grid);
		ReadVertex(name_file_vertex, vertexs);

	}


	vector<Vector3> directions;
	vector<Type> squares;

	_clock = -omp_get_wtime();
	ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, directions, squares, square_surface);
	_clock += omp_get_wtime();
	std::cout << "\n Reading time of the sphere_direction file: " << _clock << "\n";

	std::vector<cell> nodes_value;
	{
		std::vector<int> all_pairs_face;
		/*_clock = -omp_get_wtime();
		int count_unique_face = FindNeighborsPairFace(unstructured_grid, all_pairs_face);

		InitNodesValue(all_pairs_face, nodes_value);

		_clock += omp_get_wtime();*/

		FindNeighborsPairFaceAndBoundaries(unstructured_grid, all_pairs_face);

		InitNodesValue(all_pairs_face, nodes_value);

		std::unique_ptr<FILE, int(*)(FILE*)> file_id_neighbors(fopen((std::string(BASE_ADRESS) + "id_neighbors.bin").c_str(), "wb"), fclose);
		if (!file_id_neighbors) printf("id_neighbors not open\n");

		fwrite_unlocked(all_pairs_face.data(), sizeof(int), all_pairs_face.size(), file_id_neighbors.get());

		fclose(file_id_neighbors.get());
	}
	std::cout << "\n Finding time of the all_pairs_face: " << _clock << "\n";


	const int count_directions = directions.size();
	const int count_cells = unstructured_grid->GetNumberOfCells();

	InitGlobalValue(start_point_plane_coord, transform_matrix, inverse_transform_matrix, straight_face, inclined_face);


	Eigen::Matrix4d vertex_tetra;
	/* x1 x2 x3 x4
	*  y1 y2 y3 y4
	*  z1 z2 z3 z4
	*  1  1  1  1
	*/

	std::vector<Type> Illum;// (count_directions * count_cells, 0);
	std::vector<Type> Illum2;// (count_directions * count_cells, 0);


	std::vector<Normals> normals;
	ReadNormalFile(name_file_normals, normals);

	// ”пор€доченные индексы €чеек по данному направлению
	vector<int> sorted_id_cell(count_cells * count_directions);
	{
		std::unique_ptr<FILE, int(*)(FILE*)> file_id(fopen((std::string(BASE_ADRESS) + "id_defining_faces.bin").c_str(), "rb"), fclose);
		std::unique_ptr<FILE, int(*)(FILE*)>file_graph(fopen((std::string(BASE_ADRESS) + "graph" + ".bin").c_str(), "rb"), fclose);
		std::unique_ptr<FILE, int(*)(FILE*)>file_x_defining_faces(fopen((std::string(BASE_ADRESS) + "x_defining_faces" + ".bin").c_str(), "rb"), fclose);

		if (!file_id) { printf("Error file_id\n"); return 1; }
		if (!file_graph) { printf("file_graph not open\n"); return 1; }
		if (!file_x_defining_faces) { printf("file_x_defining_faces not open\n"); return 1; }

		int n;
		std::ifstream ifile(std::string(BASE_ADRESS) + "Size.txt");
		if (!ifile.is_open()) {
			printf("Error read files size\n");
			return 1;
		}
		ifile >> n;
		ifile.close();

		id_try_surface.resize(n);
		x_try_surface.resize(n);

		fread_unlocked(id_try_surface.data(), sizeof(int), n, file_id.get());
		fread_unlocked(x_try_surface.data(), sizeof(Vector3), n, file_x_defining_faces.get());
		fread_unlocked(sorted_id_cell.data(), sizeof(int), count_cells * count_directions, file_graph.get());

		fclose(file_id.get());
		fclose(file_graph.get());
		fclose(file_x_defining_faces.get());
	}


	int count = 0;
	Type norm = 0;
	Vector3 direction;

	ofstream ofile;
	ofile.open("File_with_Logs.txt");


	posX = 0;
	posX0 = 0;
	posOutC = 0;
	posOut = 0;
	posIn = 0;
	posS = 0;
	posRes = 0;
	pos_x_try = 0;

	int shift_out = 0;
	int shift_x0 = 0;
	int shift_res = 0;

	_clock = -omp_get_wtime();
	{
		std::unique_ptr<FILE, int(*)(FILE*)>file_out_id(fopen((std::string(BASE_ADRESS) + "OutId" + ".bin").c_str(), "wb"), fclose);
		std::unique_ptr<FILE, int(*)(FILE*)>file_in_id(fopen((std::string(BASE_ADRESS) + "InId" + ".bin").c_str(), "wb"), fclose);
		std::unique_ptr<FILE, int(*)(FILE*)>file_x(fopen((std::string(BASE_ADRESS) + "X" + ".bin").c_str(), "wb"), fclose);
		std::unique_ptr<FILE, int(*)(FILE*)>file_x0_local(fopen((std::string(BASE_ADRESS) + "locX0" + ".bin").c_str(), "wb"), fclose);
		std::unique_ptr<FILE, int(*)(FILE*)>file_s(fopen((std::string(BASE_ADRESS) + "S" + ".bin").c_str(), "wb"), fclose);
		std::unique_ptr<FILE, int(*)(FILE*)>file_count_out_id(fopen((std::string(BASE_ADRESS) + "CountOutId" + ".bin").c_str(), "wb"), fclose);
		std::unique_ptr<FILE, int(*)(FILE*)>file_res_bound(fopen((std::string(BASE_ADRESS) + "ResBound" + ".bin").c_str(), "wb"), fclose);

		std::unique_ptr<FILE, int(*)(FILE*)>file_shift_out(fopen((std::string(BASE_ADRESS) + "ShiftOut" + ".bin").c_str(), "wb"), fclose);
		std::unique_ptr<FILE, int(*)(FILE*)>file_shift_res(fopen((std::string(BASE_ADRESS) + "ShiftRes" + ".bin").c_str(), "wb"), fclose);
		std::unique_ptr<FILE, int(*)(FILE*)>file_shift_x0(fopen((std::string(BASE_ADRESS) + "ShiftX0" + ".bin").c_str(), "wb"), fclose);

		if (!file_out_id) printf("Outid not open\n");
		if (!file_in_id) printf("Inid not open\n");
		if (!file_x) printf("Xid not open\n");
		if (!file_x0_local) printf("X0id not open\n");
		if (!file_s) printf("Sid not open\n");
		if (!file_count_out_id) printf("file_count_out_id not open\n");
		if (!file_res_bound) printf("file_res_bound not open\n");

		if (!file_shift_out) printf("file_shift_out not open\n");
		if (!file_shift_res) printf("file_shift_res not open\n");
		if (!file_shift_x0) printf("file_shift_x0 not open\n");

		Vector3 x;
		Vector3 x0;
		int num_cell;
		const int count_cells = unstructured_grid->GetNumberOfCells();



		/*---------------------------------- далее FOR по направлени€м----------------------------------*/
		for (int num_direction = 0; num_direction < count_directions; ++num_direction)
		{
			direction = directions[num_direction];

			num_cur_direction = num_direction;
			cur_direction = direction;

			int face_state[4];  // 0=> выход€ща€ грань,  1=> вход€ща€   
		

			/*---------------------------------- далее FOR по €чейкам----------------------------------*/
			for (int h = 0; h < count_cells; ++h) {
				if (count_cells * num_direction + h > sorted_id_cell.size()) printf("err size graph\n");
				num_cell = sorted_id_cell[count_cells * num_direction + h];

				//SetVertexMatrix(num_cell, unstructured_grid, vertex_tetra);
				FindInAndOutFaces(direction, num_cell, normals, face_state);

				n_out = 0;

				for (ShortId num_out_face = 0; num_out_face < 4; ++num_out_face) {
					if (!face_state[num_out_face]) {  // выход€щие грани
						GetNodes(num_cell,grid, num_out_face, vertexs[num_cell],
							face_state, direction, normals, nodes_value,
							file_res_bound, file_s, file_x, file_x0_local, file_in_id);

						fwrite_unlocked(&num_out_face, sizeof(ShortId), 1, file_out_id.get()); posOut++;
						n_out++;
					}
				}

				fwrite_unlocked(&n_out, sizeof(ShortId), 1, file_count_out_id.get()); posOutC++;
			}
			/*---------------------------------- конец FOR по €чейкам----------------------------------*/


			fwrite_unlocked(&shift_out, sizeof(int), 1, file_shift_out.get());  shift_out = posOut;
			fwrite_unlocked(&shift_res, sizeof(int), 1, file_shift_res.get()); shift_res = posRes;
			fwrite_unlocked(&shift_x0, sizeof(int), 1, file_shift_x0.get()); shift_x0 = posX0;

			printf("End direction number #%d\n", num_direction);
		}
		/*---------------------------------- конец FOR по направлени€м----------------------------------*/

		fclose(file_res_bound.get());

		fclose(file_out_id.get());
		fclose(file_in_id.get());

		fclose(file_x.get());
		fclose(file_x0_local.get());
		fclose(file_s.get());
		fclose(file_count_out_id.get());

		fclose(file_shift_out.get());
		fclose(file_shift_res.get());
		fclose(file_shift_x0.get());

		std::ofstream ofile2(std::string(BASE_ADRESS) + "Size.txt", std::ios::app);

		ofile2 << posX << '\n';
		ofile2 << posX0 << '\n';
		ofile2 << posOutC << '\n';
		ofile2 << posOut << '\n';
		ofile2 << posIn << '\n';
		ofile2 << posS << '\n';
		ofile2 << posRes << '\n';
		/*
			X=3a
			x0
			OutC=N*M
			out=3a
			In=3a
			S=3a
			Res??
		*/

		ofile2.close();

	}


	_clock += omp_get_wtime();
	std::cout << "Time of iter: " << _clock << '\n';

	ofile << "Time of iter: " << _clock << '\n';
	ofile.close();

	return 0;
}


