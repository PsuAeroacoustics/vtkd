#include <vtkPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkCellType.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkAppendFilter.h>
#include <assert.h>

//#define _VTK_FLOAT VTK##FLOAT
//#define _VTK_DOUBLE VTK##DOUBLE

//#define _VTK(TYPE) VTK_##TYPE

namespace vtkd {

	int VTK__FLOAT() {
		return VTK_FLOAT;
	}

	int VTK__DOUBLE() {
		return VTK_DOUBLE;
	}

	int VTK__LINE() {
		return VTK_LINE;
	}

	int VTK__POLY_LINE() {
		return VTK_POLY_LINE;
	}

	int VTK__POLYGON() {
		return VTK_POLYGON;
	}

	int VTK__POLYHEDRON() {
		return VTK_POLYHEDRON;
	}

	int VTK__VOXEL() {
		return VTK_VOXEL;
	}

	class vtkPoints {
	protected:
		vtkPoints() {
			points = ::vtkPoints::New();
		}

		vtkPoints(int dataType) {
			points = ::vtkPoints::New(dataType);
		}

	public:
		::vtkPoints* points;
		void Initialize();
		bool Allocate(long size);
		vtkIdType InsertNextPoint(const double x, const double y, const double z);
		void SetPoint(long id, const double x, const double y, const double z);
		void InsertPoint(long id, const double x, const double y, const double z);
		::vtkPoints* GetPoints();
		static vtkPoints* New();
		static vtkPoints* New(int dataType);
	};


	void vtkPoints::Initialize() {
		assert(points != nullptr);
		points->Initialize();
	}

	bool vtkPoints::Allocate(long size) {
		assert(points != nullptr);
		return points->Allocate(size);
	}

	vtkIdType vtkPoints::InsertNextPoint(const double x, const double y, const double z) {
		assert(points != nullptr);
		double xyz[3] = {x, y, z};
		return points->InsertNextPoint(xyz);
	}

	void vtkPoints::SetPoint(long id, const double x, const double y, const double z) {
		assert(points != nullptr);
		points->SetPoint(id, x, y, z);
	}
	
	void vtkPoints::InsertPoint(long id, const double x, const double y, const double z) {
		assert(points != nullptr);
		points->InsertPoint(id, x, y, z);
	}

	vtkPoints* vtkPoints::New() {
		return new vtkPoints();
	}

	vtkPoints* vtkPoints::New(int dataType) {
		return new vtkPoints(dataType);
	}

	class vtkStructuredGrid {
	public:
		::vtkStructuredGrid* structured_grid;
	protected:
		vtkStructuredGrid() {
			structured_grid = ::vtkStructuredGrid::New();
		}

	public:
		void SetPoints(vtkPoints* points);
		void SetDimensions(int i, int j, int k);
		
		static vtkStructuredGrid* New();
	};

	void vtkStructuredGrid::SetPoints(vtkPoints* points) {
		assert(structured_grid != nullptr);
		structured_grid->SetPoints(points->points);
	}

	void vtkStructuredGrid::SetDimensions(int i, int j, int k) {
		structured_grid->SetDimensions(i, j, k);
	}

	vtkStructuredGrid* vtkStructuredGrid::New() {
		return new vtkStructuredGrid();
	}

	class vtkDoubleArray {
	public:
		::vtkDoubleArray* double_array;

	protected:
		vtkDoubleArray() {
			double_array = ::vtkDoubleArray::New();
		}

	public:
		void SetNumberOfComponents(int num);
		void SetNumberOfTuples(long num);
		void SetName(const char* name);
		void SetComponentName(long component, const char* name);
		void InsertTuple1(long tupleIdx, double value);
		void InsertTuple3(long tupleIdx, double val0, double val1, double val2);

		void SetTuple1(long tupleIdx, double value);
		void SetTuple3(long tupleIdx, double val0, double val1, double val2);

		static vtkDoubleArray* New();
	};

	void vtkDoubleArray::SetNumberOfComponents(int num) {
		assert(double_array != nullptr);
		double_array->SetNumberOfComponents(num);
	}

	void vtkDoubleArray::SetNumberOfTuples(long num) {
		assert(double_array != nullptr);
		double_array->SetNumberOfTuples(num);
	}

	void vtkDoubleArray::SetName(const char* name) {
		assert(double_array != nullptr);
		double_array->SetName(name);
	}

	void vtkDoubleArray::SetComponentName(long component, const char* name) {
		assert(double_array != nullptr);
		double_array->SetComponentName(component, name);
	}


	void vtkDoubleArray::InsertTuple1(long tupleIdx, double value) {
		assert(double_array != nullptr);
		double_array->InsertTuple1(tupleIdx, value);
	}

	void vtkDoubleArray::InsertTuple3(long tupleIdx, double val0, double val1, double val2) {
		assert(double_array != nullptr);
		double_array->InsertTuple3(tupleIdx, val0, val1, val2);
	}

	void vtkDoubleArray::SetTuple1(long tupleIdx, double value) {
		assert(double_array != nullptr);
		double_array->SetTuple1(tupleIdx, value);
	}

	void vtkDoubleArray::SetTuple3(long tupleIdx, double val0, double val1, double val2) {
		assert(double_array != nullptr);
		double_array->SetTuple3(tupleIdx, val0, val1, val2);
	}

	vtkDoubleArray* vtkDoubleArray::New() {
		return new vtkDoubleArray();
	}

	class vtkDataArray {
	public:
		::vtkDataArray* data_array;

		//vtkDataArray(::vtkDataArray* _data_array) : data_array(_data_array) {}
		vtkDataArray(int dataType) {
			data_array = ::vtkDataArray::CreateDataArray(dataType);
		}

		void InsertTuple1(long tupleIdx, double value);
		void InsertTuple3(long tupleIdx, double val0, double val1, double val2);

		static vtkDataArray* CreateDataArray(int dataType);
	};

	void vtkDataArray::InsertTuple1(long tupleIdx, double value) {
		assert(data_array != nullptr);
		data_array->InsertTuple1(tupleIdx, value);
	}

	void vtkDataArray::InsertTuple3(long tupleIdx, double val0, double val1, double val2) {
		assert(data_array != nullptr);
		data_array->InsertTuple3(tupleIdx, val0, val1, val2);
	}


	vtkDataArray* vtkDataArray::CreateDataArray(int dataType) {
		return new vtkDataArray(dataType);
	}

	class vtkCellData {
	public:
		::vtkCellData* cell_data;
	
		vtkCellData(::vtkCellData* _cell_data) : cell_data(_cell_data) {}

		int AddArray(vtkDoubleArray* data_array);
		int SetScalars(vtkDoubleArray* data_array);
		int SetVectors(vtkDoubleArray* data_array);
	};

	int vtkCellData::SetScalars(vtkDoubleArray* data_array) {
		assert(cell_data != nullptr);
		return cell_data->SetScalars(data_array->double_array);
	}

	int vtkCellData::SetVectors(vtkDoubleArray* data_array) {
		assert(cell_data != nullptr);
		return cell_data->SetVectors(data_array->double_array);
	}

	int vtkCellData::AddArray(vtkDoubleArray* data_array) {
		assert(cell_data != nullptr);
		return cell_data->AddArray(data_array->double_array);
	}

	class vtkPointData {
	public:
		::vtkPointData* point_data;
	
		vtkPointData(::vtkPointData* _point_data) : point_data(_point_data) {}

		int AddArray(vtkDoubleArray* data_array);
		int SetScalars(vtkDoubleArray* data_array);
		int SetVectors(vtkDoubleArray* data_array);
	};

	int vtkPointData::SetScalars(vtkDoubleArray* data_array) {
		assert(point_data != nullptr);
		return point_data->SetScalars(data_array->double_array);
	}

	int vtkPointData::SetVectors(vtkDoubleArray* data_array) {
		assert(point_data != nullptr);
		return point_data->SetVectors(data_array->double_array);
	}

	int vtkPointData::AddArray(vtkDoubleArray* data_array) {
		assert(point_data != nullptr);
		return point_data->AddArray(data_array->double_array);
	}

	class vtkUnstructuredGrid {
	public:
		::vtkUnstructuredGrid* grid;
		vtkUnstructuredGrid(::vtkUnstructuredGrid* _grid) : grid(_grid) {}
	protected:
		vtkUnstructuredGrid() {
			grid = ::vtkUnstructuredGrid::New();
		}

	public:
		void SetPoints(vtkPoints* points);
		long InsertNextCell(int type, long npts, long* ptIds);
		vtkIdType FindPoint(double x, double y, double z);
		vtkPointData* GetPointData();
		vtkCellData* GetCellData();
		void Allocate(long numCells);
		static vtkUnstructuredGrid* New();
	};

	vtkUnstructuredGrid* vtkUnstructuredGrid::New() {
		return new vtkUnstructuredGrid();
	}

	void vtkUnstructuredGrid::SetPoints(vtkPoints* points) {
		assert(grid != nullptr);
		grid->SetPoints(points->points);
	}

	long vtkUnstructuredGrid::InsertNextCell(int type, long npts, long* ptIds) {
		assert(grid != nullptr);
		return grid->InsertNextCell(type, static_cast<vtkIdType>(npts), reinterpret_cast<vtkIdType*>(ptIds));
	}

	vtkIdType vtkUnstructuredGrid::FindPoint(double x, double y, double z) {
		assert(grid != nullptr);
		return grid->FindPoint(x, y, z);
	}

	vtkCellData* vtkUnstructuredGrid::GetCellData() {
		assert(grid != nullptr);
		auto cd = grid->GetCellData();
		assert(cd != nullptr);
		return new vtkCellData(cd);
	}

	vtkPointData* vtkUnstructuredGrid::GetPointData() {
		assert(grid != nullptr);
		auto pd = grid->GetPointData();
		assert(pd != nullptr);
		return new vtkPointData(pd);
	}

	void vtkUnstructuredGrid::Allocate(long numCells) {
		assert(grid != nullptr);
		grid->Allocate(numCells);
	}

	class vtkXMLStructuredGridWriter {
	private:
		::vtkXMLStructuredGridWriter* writer;
	protected:
		vtkXMLStructuredGridWriter() {
			writer = ::vtkXMLStructuredGridWriter::New();
		}

	public:
		void SetFileName(const char* filename);
		void SetInputData(vtkStructuredGrid* data);
		int Write();

		static vtkXMLStructuredGridWriter* New();
	};

	void vtkXMLStructuredGridWriter::SetFileName(const char* filename) {
		writer->SetFileName(filename);
	}

	void vtkXMLStructuredGridWriter::SetInputData(vtkStructuredGrid* data) {
		writer->SetInputData(data->structured_grid);
	}

	int vtkXMLStructuredGridWriter::Write() {
		return writer->Write();
	}
	
	vtkXMLStructuredGridWriter* vtkXMLStructuredGridWriter::New() {
		return new vtkXMLStructuredGridWriter();
	}



	class vtkXMLUnstructuredGridWriter {
	private:
		::vtkXMLUnstructuredGridWriter* writer;
	protected:
		vtkXMLUnstructuredGridWriter() {
			writer = ::vtkXMLUnstructuredGridWriter::New();
		}

	public:
		void SetFileName(const char* filename);
		void SetInputData(vtkUnstructuredGrid* data);
		void SetInputData(int i, vtkUnstructuredGrid* data);
		void SetDataModeToAscii();
		int Write();

		static vtkXMLUnstructuredGridWriter* New();
	};

	void vtkXMLUnstructuredGridWriter::SetFileName(const char* filename) {
		assert(writer != nullptr);
		writer->SetFileName(filename);
	}

	void vtkXMLUnstructuredGridWriter::SetDataModeToAscii() {
		assert(writer != nullptr);
		writer->SetDataModeToAscii();
	}

	void vtkXMLUnstructuredGridWriter::SetInputData(vtkUnstructuredGrid* data) {
		assert(writer != nullptr);
		writer->SetInputData(data->grid);
		//writer->SetInputData(data);
	}

	void vtkXMLUnstructuredGridWriter::SetInputData(int i, vtkUnstructuredGrid* data) {
		assert(writer != nullptr);
		writer->SetInputData(i, data->grid);
		//writer->SetInputData(data);
	}

	int vtkXMLUnstructuredGridWriter::Write() {
		assert(writer != nullptr);
		return writer->Write();
	}
	
	vtkXMLUnstructuredGridWriter* vtkXMLUnstructuredGridWriter::New() {
		return new vtkXMLUnstructuredGridWriter();
	}



	class vtkXMLPUnstructuredGridWriter {
	private:
		::vtkXMLPUnstructuredGridWriter* writer;
	protected:
		vtkXMLPUnstructuredGridWriter() {
			writer = ::vtkXMLPUnstructuredGridWriter::New();
		}

	public:
		void SetFileName(const char* filename);
		void SetInputData(vtkUnstructuredGrid* data);
		int Write();
		void SetNumberOfPieces(int num);
		void SetStartPiece(int start);
		void SetEndPiece(int end);

		static vtkXMLPUnstructuredGridWriter* New();
	};

	void vtkXMLPUnstructuredGridWriter::SetFileName(const char* filename) {
		assert(writer != nullptr);
		writer->SetFileName(filename);
	}

	void vtkXMLPUnstructuredGridWriter::SetInputData(vtkUnstructuredGrid* data) {
		assert(writer != nullptr);
		writer->SetInputData(data->grid);
	}

	void vtkXMLPUnstructuredGridWriter::SetNumberOfPieces(int num) {
		assert(writer != nullptr);
		writer->SetNumberOfPieces(num);
	}

	int vtkXMLPUnstructuredGridWriter::Write() {
		assert(writer != nullptr);
		return writer->Write();
	}
	
	void vtkXMLPUnstructuredGridWriter::SetStartPiece(int start) {
		assert(writer != nullptr);
		writer->SetStartPiece(start);
	}

	void vtkXMLPUnstructuredGridWriter::SetEndPiece(int end) {
		assert(writer != nullptr);
		writer->SetEndPiece(end);
	}


	vtkXMLPUnstructuredGridWriter* vtkXMLPUnstructuredGridWriter::New() {
		return new vtkXMLPUnstructuredGridWriter();
	}

	class vtkAppendFilter {
	protected:
		vtkAppendFilter() {
			append_filter = ::vtkAppendFilter::New();
		}
	public:
		::vtkAppendFilter* append_filter;

		void Update();
		void AddInputData(vtkUnstructuredGrid* data);
		void AddInputData(int i, vtkUnstructuredGrid* data);

		void SetInputData(vtkUnstructuredGrid* data);
		void SetInputData(int i, vtkUnstructuredGrid* data);

		void RemoveAllInputs();
		vtkUnstructuredGrid* GetOutput();
		static vtkAppendFilter* New();
	};

	void vtkAppendFilter::Update() {
		assert(append_filter != nullptr);
		append_filter->Update();
	}

	void vtkAppendFilter::AddInputData(vtkUnstructuredGrid* data) {
		assert(append_filter != nullptr);
		append_filter->AddInputData(data->grid);
	}

	void vtkAppendFilter::AddInputData(int i, vtkUnstructuredGrid* data) {
		assert(append_filter != nullptr);
		append_filter->AddInputData(i, data->grid);
	}

	void vtkAppendFilter::SetInputData(vtkUnstructuredGrid* data) {
		assert(append_filter != nullptr);
		append_filter->SetInputData(data->grid);
	}

	void vtkAppendFilter::SetInputData(int i, vtkUnstructuredGrid* data) {
		assert(append_filter != nullptr);
		append_filter->SetInputData(i, data->grid);
	}

	void vtkAppendFilter::RemoveAllInputs() {
		assert(append_filter != nullptr);
		append_filter->RemoveAllInputs();

	}

	vtkUnstructuredGrid* vtkAppendFilter::GetOutput() {
		assert(append_filter != nullptr);
		return new vtkUnstructuredGrid(append_filter->GetOutput());
	}

	vtkAppendFilter* vtkAppendFilter::New() {
		return new vtkAppendFilter();
	}
}
