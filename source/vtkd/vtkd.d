module vtkd.vtkd;

extern(C++, vtkd) {
	
	alias vtkIdType = long;

	int VTK__FLOAT();
	int VTK__DOUBLE();
	int VTK__LINE();
	int VTK__POLY_LINE();
	int VTK__POLYGON();
	int VTK__POLYHEDRON();
	int VTK__VOXEL();

	extern(C++, class) struct vtkAppendFilter {
		@nogc void Update();
		@nogc void AddInputData(vtkUnstructuredGrid* data);
		@nogc void AddInputData(int i, vtkUnstructuredGrid* data);
		@nogc void SetInputData(vtkUnstructuredGrid* data);
		@nogc void SetInputData(int i, vtkUnstructuredGrid* data);
		@nogc void RemoveAllInputs();
		@nogc vtkUnstructuredGrid* GetOutput();
		@nogc static vtkAppendFilter* New();
	}

	extern(C++, class) struct vtkPoints {
		@nogc void Initialize();
		@nogc bool Allocate(vtkIdType size);
		@nogc vtkIdType InsertNextPoint(const double x, const double y, const double z);
		@nogc void SetPoint(vtkIdType id, const double x, const double y, const double z);
		@nogc void InsertPoint(vtkIdType id, const double x, const double y, const double z);

		@nogc static vtkPoints* New();
		@nogc static vtkPoints* New(int dataType);
	}

	extern(C++, class) struct vtkStructuredGrid {
		@nogc void SetPoints(vtkPoints* points);
		@nogc void SetDimensions(int i, int j, int k);

		@nogc static vtkStructuredGrid* New();
	}


	extern(C++, class) struct vtkUnstructuredGrid {
		@nogc void SetPoints(vtkPoints* points);
		@nogc void Allocate(vtkIdType numCells);
		@nogc vtkIdType InsertNextCell(int type, vtkIdType npts, vtkIdType* ptIds);
		@nogc vtkIdType FindPoint(double x, double y, double z);
		@nogc vtkCellData* GetCellData();
		@nogc vtkPointData* GetPointData();
		@nogc static vtkUnstructuredGrid* New();
	}


	extern(C++, class) struct vtkXMLStructuredGridWriter {
		@nogc void SetFileName(const char* filename);
		@nogc void SetInputData(vtkStructuredGrid* data);
		@nogc int Write();

		@nogc static vtkXMLStructuredGridWriter* New();
	}
	
	extern(C++, class) struct vtkXMLUnstructuredGridWriter {

		@nogc void SetFileName(const char* filename);
		@nogc void SetInputData(vtkUnstructuredGrid* data);
		@nogc void SetInputData(int i, vtkUnstructuredGrid* data);
		@nogc void SetDataModeToAscii();
		@nogc int Write();

		@nogc static vtkXMLUnstructuredGridWriter* New();
	}

	extern(C++, class) struct vtkXMLPUnstructuredGridWriter {
		@nogc void SetFileName(const char* filename);
		@nogc void SetInputData(vtkUnstructuredGrid* data);
		@nogc int Write();
		@nogc void SetNumberOfPieces(int num);
		@nogc void SetStartPiece(int start);
		@nogc void SetEndPiece(int end);

		@nogc static vtkXMLPUnstructuredGridWriter* New();
	}

	extern(C++, class) struct vtkDataArray {
		@nogc void InsertTuple1(vtkIdType tupleIdx, double value);
		@nogc void InsertTuple3(vtkIdType tupleIdx, double val0, double val1, double val2);

		@nogc static vtkDataArray* CreateDataArray(int dataType);
	}

	extern(C++, class) struct vtkDoubleArray {
		@nogc void SetNumberOfComponents(int num);
		@nogc void SetNumberOfTuples(vtkIdType num);
		@nogc void SetName(const char* name);
		@nogc void SetComponentName(long component, const char* name);
		@nogc void InsertTuple1(vtkIdType tupleIdx, double value);
		@nogc void InsertTuple3(vtkIdType tupleIdx, double val0, double val1, double val2);
		@nogc void SetTuple1(vtkIdType tupleIdx, double value);
		@nogc void SetTuple3(vtkIdType tupleIdx, double val0, double val1, double val2);
		@nogc static vtkDoubleArray* New();
	}

	extern(C++, class) struct vtkCellData {
		@nogc int AddArray(vtkDoubleArray* data_array);
		@nogc int SetScalars(vtkDoubleArray* data_array);
		@nogc int SetVectors(vtkDoubleArray* data_array);
	}

	extern(C++, class) struct vtkPointData {
		@nogc int AddArray(vtkDoubleArray* data_array);
		@nogc int SetScalars(vtkDoubleArray* data_array);
		@nogc int SetVectors(vtkDoubleArray* data_array);
	}
}


import albm.immersedboundary;

auto save_to_vtu(SimConfig, VelField, DenField)(auto ref SimConfig config, auto ref VelField velocity_field, auto ref DenField density_field, string filename) {
	
	vtkIdType[SimConfig.PosVector] node_id_map;

	auto unstructured_grid = vtkUnstructuredGrid.New;
	auto points = vtkPoints.New(VTK__DOUBLE);
	auto grid_writer = vtkXMLUnstructuredGridWriter.New;

	points.Allocate(config.lattice.z_nodes*config.lattice.x_nodes*config.lattice.y_nodes);

	foreach(z_block, y_list; config.immersed_boundary.regular_blocks) {
		foreach(y_block, x_list; y_list) {
			foreach(x_block; x_list) {
				foreach(chunk_idx; 0..SimConfig.BlockType.total_chunks) {
					auto pos = SimConfig.BlockType.get_pos(chunk_idx);
					auto y = pos[1] + SimConfig.BlockType.y_nodes*y_block;
					auto z = pos[2] + SimConfig.BlockType.z_nodes*z_block;

					foreach(node_idx; 0..SimConfig.chunk_size) {
						auto x = pos[0]*SimConfig.chunk_size + SimConfig.BlockType.x_nodes*x_block + node_idx;

						auto node_pos = SimConfig.PosVector(x, y, z);
						node_id_map[node_pos] = points.InsertNextPoint(x, y, z);//x + y*config.lattice.x_nodes + z*config.lattice.x_nodes*config.lattice.y_nodes;
					}
				}
			}
		}
	}

	foreach(z_block, y_list; config.immersed_boundary.immersed_blocks) {
		foreach(y_block, x_list; y_list) {
			foreach(x_block; x_list) {

				auto addr = LBAddress(x_block, y_block, z_block);
				foreach(chunk_idx, ref link_block; config.immersed_boundary.link_blocks[addr]) {
				
					auto pos = SimConfig.BlockType.get_pos(chunk_idx);
					auto y = pos[1] + SimConfig.BlockType.y_nodes*y_block;
					auto z = pos[2] + SimConfig.BlockType.z_nodes*z_block;

					foreach(node_idx; 0..SimConfig.chunk_size) {
						if(link_block.node_type[node_idx] != NodeType.solid) {
							auto x = pos[0]*SimConfig.chunk_size + SimConfig.BlockType.x_nodes*x_block + node_idx;
							auto node_pos = SimConfig.PosVector(x, y, z);
							//node_id_map[node_pos] = x + y*config.lattice.x_nodes + z*config.lattice.x_nodes*config.lattice.y_nodes;
							node_id_map[node_pos] = points.InsertNextPoint(x, y, z);
						}
					}
				}
			}
		}
	}

	foreach(z; 0..config.lattice.z_nodes) {
		foreach(y; 0..config.lattice.y_nodes) {
			foreach(x; 0..config.lattice.x_nodes) {
				auto pos1 = SimConfig.PosVector(x, y, z);
				auto pos2 = SimConfig.PosVector(x + 1, y, z);
				auto pos3 = SimConfig.PosVector(x, y + 1, z);
				auto pos4 = SimConfig.PosVector(x + 1, y + 1, z);
				
				auto pos5 = SimConfig.PosVector(x, y, z + 1);
				auto pos6 = SimConfig.PosVector(x + 1, y, z + 1);
				auto pos7 = SimConfig.PosVector(x, y + 1, z + 1);
				auto pos8 = SimConfig.PosVector(x + 1, y + 1, z + 1);
				
				if((pos1 in node_id_map) && (pos2 in node_id_map) &&
					(pos3 in node_id_map) && (pos4 in node_id_map) &&
					(pos5 in node_id_map) && (pos6 in node_id_map) &&
					(pos7 in node_id_map) && (pos8 in node_id_map)

				) {
					vtkIdType[8] ids = [
						node_id_map[pos1],
						node_id_map[pos2],
						node_id_map[pos3],
						node_id_map[pos4],
						node_id_map[pos5],
						node_id_map[pos6],
						node_id_map[pos7],
						node_id_map[pos8]
					];

					unstructured_grid.InsertNextCell(VTK__VOXEL, ids.length, ids.ptr);
				}
			}
		}
	}

	auto density_array = vtkDoubleArray.New();
	density_array.SetNumberOfComponents(1);
	density_array.SetNumberOfTuples(config.lattice.total_nodes);
	density_array.SetName("density");

	auto velocity_array = vtkDoubleArray.New();
	velocity_array.SetNumberOfComponents(3);
	velocity_array.SetNumberOfTuples(config.lattice.total_nodes);
	velocity_array.SetName("velocity");
	velocity_array.SetComponentName(0, "u");
	velocity_array.SetComponentName(1, "v");
	velocity_array.SetComponentName(2, "w");

	foreach(pos, id; node_id_map) {
		//auto id = idx - mpi_manager.offset[2]*mpi_manager.lattice.x_nodes*mpi_manager.lattice.y_nodes;
		density_array.SetTuple1(id, density_field[pos[0], pos[1], pos[2]]);
		auto vel = velocity_field[pos[0], pos[1], pos[2]];
		velocity_array.SetTuple3(id, vel[0], vel[1], vel[2]);
	}
	
	unstructured_grid.SetPoints(points);

	auto point_data = unstructured_grid.GetPointData;
	point_data.SetScalars(density_array);
	point_data.SetVectors(velocity_array);

	import std.string : toStringz;
	//auto grid_filename = "domain.vtu";
	//grid_writer.SetDataModeToAscii;
	grid_writer.SetFileName(filename.toStringz);
	grid_writer.SetInputData(unstructured_grid);
	grid_writer.Write;
}

auto save_to_pvtu(SimConfig, MPIManager, VelField, DenField)(auto ref SimConfig config, auto ref MPIManager mpi_manager, auto ref VelField velocity_field, auto ref DenField density_field, string fileprefix) {
	
	import mpid;
	
	import numd.linearalgebra.matrix : Vector;
	
	import std.algorithm : sum;
	import std.conv : to;

	alias Vec3 = Vector!(SimConfig.StencilType.dimensions, SimConfig.ComputeType);

	alias T = SimConfig.ComputeType;

	vtkIdType[SimConfig.PosVector] node_id_map;

	auto unstructured_grid = vtkUnstructuredGrid.New;
	auto points = vtkPoints.New(VTK__DOUBLE);
	auto grid_writer = vtkXMLUnstructuredGridWriter.New;

	points.Allocate((mpi_manager.lattice.z_nodes + 1)*mpi_manager.lattice.x_nodes*mpi_manager.lattice.y_nodes);

	auto z_start_block = mpi_manager.z_blocks_per_proc[0..mpi_manager.global_comm.rank].sum;
	auto z_end_block = z_start_block + mpi_manager.z_blocks_per_proc[mpi_manager.global_comm.rank] + 1;

	if(z_end_block > config.lattice.z_blocks) {
		z_end_block = config.lattice.z_blocks;
	}

	foreach(z_block, y_list; config.immersed_boundary.regular_blocks[z_start_block..z_end_block]) {
		z_block += z_start_block;
		foreach(y_block, x_list; y_list) {
			foreach(x_block; x_list) {
				foreach(chunk_idx; 0..SimConfig.BlockType.total_chunks) {
					auto pos = SimConfig.BlockType.get_pos(chunk_idx);
					auto y = pos[1] + SimConfig.BlockType.y_nodes*y_block;
					auto z = pos[2] + SimConfig.BlockType.z_nodes*z_block;

					foreach(node_idx; 0..SimConfig.chunk_size) {
						auto x = pos[0]*SimConfig.chunk_size + SimConfig.BlockType.x_nodes*x_block + node_idx;

						auto node_pos = SimConfig.PosVector(x, y, z);
						node_id_map[node_pos] = x + y*config.lattice.x_nodes + z*config.lattice.x_nodes*config.lattice.y_nodes;
					}
				}
			}
		}
	}

	foreach(z_block, y_list; config.immersed_boundary.immersed_blocks[z_start_block..z_end_block]) {
		z_block += z_start_block;
		foreach(y_block, x_list; y_list) {
			foreach(x_block; x_list) {

				auto addr = LBAddress(x_block, y_block, z_block);
				foreach(chunk_idx, ref link_block; config.immersed_boundary.link_blocks[addr]) {
				
					auto pos = SimConfig.BlockType.get_pos(chunk_idx);
					auto y = pos[1] + SimConfig.BlockType.y_nodes*y_block;
					auto z = pos[2] + SimConfig.BlockType.z_nodes*z_block;

					foreach(node_idx; 0..SimConfig.chunk_size) {
						if(link_block.node_type[node_idx] != NodeType.solid) {
							auto x = pos[0]*SimConfig.chunk_size + SimConfig.BlockType.x_nodes*x_block + node_idx;
							auto node_pos = SimConfig.PosVector(x, y, z);
							node_id_map[node_pos] = x + y*config.lattice.x_nodes + z*config.lattice.x_nodes*config.lattice.y_nodes;
						}
					}
				}
			}
		}
	}

	foreach(z; mpi_manager.offset[2]..mpi_manager.lattice.z_nodes + mpi_manager.offset[2]) {
		foreach(y; mpi_manager.offset[1]..mpi_manager.lattice.y_nodes + mpi_manager.offset[1]) {
			foreach(x; mpi_manager.offset[0]..mpi_manager.lattice.x_nodes + mpi_manager.offset[0]) {
				auto pos1 = SimConfig.PosVector(x, y, z);
				auto pos2 = SimConfig.PosVector(x + 1, y, z);
				auto pos3 = SimConfig.PosVector(x, y + 1, z);
				auto pos4 = SimConfig.PosVector(x + 1, y + 1, z);
				
				auto pos5 = SimConfig.PosVector(x, y, z + 1);
				auto pos6 = SimConfig.PosVector(x + 1, y, z + 1);
				auto pos7 = SimConfig.PosVector(x, y + 1, z + 1);
				auto pos8 = SimConfig.PosVector(x + 1, y + 1, z + 1);
				
				if((pos1 in node_id_map) && (pos2 in node_id_map) &&
					(pos3 in node_id_map) && (pos4 in node_id_map) &&
					(pos5 in node_id_map) && (pos6 in node_id_map) &&
					(pos7 in node_id_map) && (pos8 in node_id_map)

				) {
					vtkIdType[8] ids = [
						node_id_map[pos1] - mpi_manager.offset[2]*mpi_manager.lattice.x_nodes*mpi_manager.lattice.y_nodes,
						node_id_map[pos2] - mpi_manager.offset[2]*mpi_manager.lattice.x_nodes*mpi_manager.lattice.y_nodes,
						node_id_map[pos3] - mpi_manager.offset[2]*mpi_manager.lattice.x_nodes*mpi_manager.lattice.y_nodes,
						node_id_map[pos4] - mpi_manager.offset[2]*mpi_manager.lattice.x_nodes*mpi_manager.lattice.y_nodes,
						node_id_map[pos5] - mpi_manager.offset[2]*mpi_manager.lattice.x_nodes*mpi_manager.lattice.y_nodes,
						node_id_map[pos6] - mpi_manager.offset[2]*mpi_manager.lattice.x_nodes*mpi_manager.lattice.y_nodes,
						node_id_map[pos7] - mpi_manager.offset[2]*mpi_manager.lattice.x_nodes*mpi_manager.lattice.y_nodes,
						node_id_map[pos8] - mpi_manager.offset[2]*mpi_manager.lattice.x_nodes*mpi_manager.lattice.y_nodes
					];

					unstructured_grid.InsertNextCell(VTK__VOXEL, ids.length, ids.ptr);
					points.InsertPoint(ids[0], pos1[0], pos1[1], pos1[2]);
					points.InsertPoint(ids[1], pos2[0], pos2[1], pos2[2]);
					points.InsertPoint(ids[2], pos3[0], pos3[1], pos3[2]);
					points.InsertPoint(ids[3], pos4[0], pos4[1], pos4[2]);
					points.InsertPoint(ids[4], pos5[0], pos5[1], pos5[2]);
					points.InsertPoint(ids[5], pos6[0], pos6[1], pos6[2]);
					points.InsertPoint(ids[6], pos7[0], pos7[1], pos7[2]);
					points.InsertPoint(ids[7], pos8[0], pos8[1], pos8[2]);
				}
			}
		}
	}

	auto density_array = vtkDoubleArray.New();
	density_array.SetNumberOfComponents(1);
	density_array.SetNumberOfTuples(mpi_manager.lattice.total_nodes + mpi_manager.lattice.x_nodes*mpi_manager.lattice.y_nodes);
	density_array.SetName("density");

	auto velocity_array = vtkDoubleArray.New();
	velocity_array.SetNumberOfComponents(3);
	velocity_array.SetNumberOfTuples(mpi_manager.lattice.total_nodes + mpi_manager.lattice.x_nodes*mpi_manager.lattice.y_nodes);
	velocity_array.SetName("velocity");
	velocity_array.SetComponentName(0, "u");
	velocity_array.SetComponentName(1, "v");
	velocity_array.SetComponentName(2, "w");

	T[] send_density = new T[mpi_manager.lattice.x_nodes * mpi_manager.lattice.y_nodes];
	Vec3[] send_velocity = new Vec3[mpi_manager.lattice.x_nodes * mpi_manager.lattice.y_nodes];
	T[] recv_density = new T[mpi_manager.lattice.x_nodes * mpi_manager.lattice.y_nodes];
	Vec3[] recv_velocity = new Vec3[mpi_manager.lattice.x_nodes * mpi_manager.lattice.y_nodes];

	foreach(y; mpi_manager.offset[1]..mpi_manager.lattice.y_nodes + mpi_manager.offset[1]) {
		foreach(x; mpi_manager.offset[0]..mpi_manager.lattice.x_nodes + mpi_manager.offset[0]) {
			send_density[x + y*mpi_manager.lattice.y_nodes] = density_field[x, y, 0];
			send_velocity[x + y*mpi_manager.lattice.y_nodes] = velocity_field[x, y, 0];
		}
	}

	if(mpi_manager.global_comm.rank > 0) {
		mpi_manager.global_comm.send_array(send_density, mpi_manager.global_comm.rank - 1, 0);
	}

	if(mpi_manager.global_comm.rank < mpi_manager.global_comm.size - 1) {
		recv_density = mpi_manager.global_comm.recv_array!T(mpi_manager.global_comm.rank + 1, 0);

		auto idx = 0;
		foreach(y; mpi_manager.offset[1]..mpi_manager.lattice.y_nodes + mpi_manager.offset[1]) {
			foreach(x; mpi_manager.offset[0]..mpi_manager.lattice.x_nodes + mpi_manager.offset[0]) {
				vtkIdType id = x + y*config.lattice.x_nodes + (mpi_manager.lattice.z_nodes)*config.lattice.x_nodes*config.lattice.y_nodes;
				density_array.SetTuple1(id, recv_density[idx]);
				idx++;
			}
		}
	}

	if(mpi_manager.global_comm.rank > 0) {
		mpi_manager.global_comm.send_array(send_velocity, mpi_manager.global_comm.rank - 1, 0);
	}

	if(mpi_manager.global_comm.rank < mpi_manager.global_comm.size - 1) {
		recv_velocity = mpi_manager.global_comm.recv_array!Vec3(mpi_manager.global_comm.rank + 1, 0);

		auto idx = 0;
		foreach(y; mpi_manager.offset[1]..mpi_manager.lattice.y_nodes + mpi_manager.offset[1]) {
			foreach(x; mpi_manager.offset[0]..mpi_manager.lattice.x_nodes + mpi_manager.offset[0]) {
				vtkIdType id = x + y*config.lattice.x_nodes + mpi_manager.lattice.z_nodes*config.lattice.x_nodes*config.lattice.y_nodes;
				velocity_array.SetTuple3(id, recv_velocity[idx][0], recv_velocity[idx][1], recv_velocity[idx][2]);
				idx++;
			}
		}
	}

	foreach(pos, idx; node_id_map) {
		auto local_pos = pos - mpi_manager.offset;
		if((local_pos[0] < mpi_manager.lattice.x_nodes) && (local_pos[1] < mpi_manager.lattice.y_nodes) && (local_pos[2] < mpi_manager.lattice.z_nodes)) {
			auto id = idx - mpi_manager.offset[2]*mpi_manager.lattice.x_nodes*mpi_manager.lattice.y_nodes;
			density_array.SetTuple1(id, density_field[local_pos[0], local_pos[1], local_pos[2]]);
			auto vel = velocity_field[local_pos[0], local_pos[1], local_pos[2]];
			velocity_array.SetTuple3(id, vel[0], vel[1], vel[2]);
		}
	}
	
	unstructured_grid.SetPoints(points);

	auto point_data = unstructured_grid.GetPointData;
	point_data.SetScalars(density_array);
	point_data.SetVectors(velocity_array);

	if(mpi_manager.global_comm.rank == 0) {
		import dxml.writer;

		import std.array : appender;
		import std.string : toStringz;
		
		// Write our own pvtu file because libvtk's is stupid.
		// Basically vtkPUnstructuredGrid writer wants to 
		// write both the pvtu file (small) and the vtu files
		// for each proc, but it does it with its local data for
		// whatever reason. This basically means if we use it
		// we would first write out the wrong data, and then have
		// every proc write out the correct data. Also for whatever
		// reason wont but the velocity component names in the pvtu
		// file, which we do.
		auto writer = xmlWriter(appender!string);

		writer.openStartTag("VTKFile");
		writer.writeAttr("type", "PUnstructuredGrid");
		writer.writeAttr("version", "0.1");
		writer.writeAttr("byte_order", "LittleEndian");
		writer.writeAttr("header_type", "UInt32");
		writer.writeAttr("compressor", "vtkZLibDataCompressor");
		writer.closeStartTag;

		writer.openStartTag("PUnstructuredGrid");
		writer.writeAttr("GhostLevel", "0");
		writer.closeStartTag;

		writer.openStartTag("PPointData");
		writer.writeAttr("Scalars", "density");
		writer.writeAttr("Vectors", "velocity");
		writer.closeStartTag;

		writer.openStartTag("PDataArray");
		writer.writeAttr("type", "Float64");
		writer.writeAttr("Name", "density");
		writer.closeStartTag(EmptyTag.yes);

		writer.openStartTag("PDataArray");
		writer.writeAttr("type", "Float64");
		writer.writeAttr("Name", "velocity");
		writer.writeAttr("NumberOfComponents", "3");
		writer.writeAttr("ComponentName0", "u");
		writer.writeAttr("ComponentName1", "v");
		writer.writeAttr("ComponentName2", "w");
		writer.closeStartTag(EmptyTag.yes);

		writer.writeEndTag("PPointData");

		writer.writeStartTag("PPoints");

		writer.openStartTag("PDataArray");
		writer.writeAttr("type", "Float64");
		writer.writeAttr("Name", "Points");
		writer.writeAttr("NumberOfComponents", "3");
		writer.closeStartTag(EmptyTag.yes);

		writer.writeEndTag("PPoints");

		foreach(peice; 0..mpi_manager.global_comm.size) {
			writer.openStartTag("Piece");
			auto grid_filename = fileprefix~"_"~peice.to!string~".vtu";
			writer.writeAttr("Source", grid_filename);
			writer.closeStartTag(EmptyTag.yes);
		}

		writer.writeEndTag("PUnstructuredGrid");
		writer.writeEndTag("VTKFile");

		import std.file : write;
		auto output_filename = fileprefix~".pvtu";
		write(output_filename, "<?xml version=\"1.0\"?>\n"~writer.output.data);
	}

	mpi_manager.global_comm.barrier;

	import std.string : toStringz;
	auto grid_filename = fileprefix~"_"~mpi_manager.global_comm.rank.to!string~".vtu";
	grid_writer.SetFileName(grid_filename.toStringz);
	grid_writer.SetInputData(unstructured_grid);
	grid_writer.Write;
}
