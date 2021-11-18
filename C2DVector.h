#ifndef C2DVECTOR_H_HEADER_INCLUDED_BECEF599
#define C2DVECTOR_H_HEADER_INCLUDED_BECEF599
#include <vector>
using namespace std;

template <class T>
class C2DVector
{
public:
	C2DVector():m_dimRow(0), m_dimCol(0){;}
	C2DVector(int nRow, int nCol) {
		m_dimRow = nRow;
		m_dimCol = nCol;
		for (int i=0; i < nRow; i++){
			vector<T> x(nCol);
			int y = x.size();
			m_2DVector.push_back(x);
		}
	}
	void SetAt(int nRow, int nCol, const T& value) throw(std::out_of_range) {
		if(nRow >= m_dimRow || nCol >= m_dimCol)
			throw out_of_range("Array out of bound");
		else
			m_2DVector[nRow][nCol] = value;
	}
	T GetAt(int nRow, int nCol) {
		if(nRow >= m_dimRow || nCol >= m_dimCol)
			throw out_of_range("Array out of bound");
		else
			return m_2DVector[nRow][nCol];
	}
	void GrowRow(int newSize) {
		if (newSize <= m_dimRow)
			return;
		m_dimRow = newSize;
		for(int i = 0 ; i < newSize - m_dimCol; i++)   {
			vector<int> x(m_dimRow);
			m_2DVector.push_back(x);
		}
	}
	void GrowCol(int newSize) {
		if(newSize <= m_dimCol)
			return;
		m_dimCol = newSize;
		for (int i=0; i <m_dimRow; i++)
			m_2DVector[i].resize(newSize);
	}
	vector<T>& operator[](int x)    {
		return m_2DVector[x];
	}
private:
	vector< vector <T> > m_2DVector;
	unsigned int m_dimRow;
	unsigned int m_dimCol;
};


#endif 
