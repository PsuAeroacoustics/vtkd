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
		@nogc void Delete();
		@nogc static vtkAppendFilter* New();
	}

	extern(C++, class) struct vtkPoints {
		@nogc void Delete();
		@nogc void Initialize();
		@nogc bool Allocate(vtkIdType size);
		@nogc vtkIdType InsertNextPoint(const double x, const double y, const double z);
		@nogc void SetPoint(vtkIdType id, const double x, const double y, const double z);
		@nogc void InsertPoint(vtkIdType id, const double x, const double y, const double z);

		@nogc static vtkPoints* New();
		@nogc static vtkPoints* New(int dataType);
	}

	extern(C++, class) struct vtkStructuredGrid {
		@nogc void Delete();
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
		@nogc void Delete();
		@nogc static vtkUnstructuredGrid* New();
	}


	extern(C++, class) struct vtkXMLStructuredGridWriter {
		@nogc void SetFileName(const char* filename);
		@nogc void SetInputData(vtkStructuredGrid* data);
		@nogc int Write();
		@nogc void Delete();

		@nogc static vtkXMLStructuredGridWriter* New();
	}
	
	extern(C++, class) struct vtkXMLUnstructuredGridWriter {
		@nogc void Delete();
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
		@nogc void Delete();
		@nogc void InsertTuple1(vtkIdType tupleIdx, double value);
		@nogc void InsertTuple3(vtkIdType tupleIdx, double val0, double val1, double val2);

		@nogc static vtkDataArray* CreateDataArray(int dataType);
	}

	extern(C++, class) struct vtkDoubleArray {
		@nogc void Delete();
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
