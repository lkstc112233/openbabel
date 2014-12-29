// TemplateHDVector.hpp
// The high degree vector.
// Min val find function included!
#pragma once

#include <vector>
template<class T>
class MultiVector
{
private:
	union
	{
		MultiVector<T>**	lowerVector;
		T*			basicVector;
	} value;
	int m_degrees;
	int m_size;
public:
	MultiVector(int degrees,int size)
		: m_degrees(degrees)
		, m_size(size)
	{
		if (1==degrees)
			value.basicVector=new T[m_size];
		else
		{
			value.lowerVector=new MultiVector<T>*[m_size];
			for (int i=0;i<m_size;++i)
				value.lowerVector[i]=new MultiVector<T>(degrees-1,m_size);
		}
	}
	~MultiVector()
	{
		if (1==m_degrees)
			delete [] value.basicVector;
		else
		{
			for (int i=0;i<m_size;++i)
				delete value.lowerVector[i];
			delete[] value.lowerVector;
		}
	}
	void fill(const T &newValue)
	{
		if (1==m_degrees)
			for (int i=0;i<m_size;++i)
				value.basicVector[i]=newValue;
		else
			for (int i=0;i<m_size;++i)
				value.lowerVector[i]->fill(newValue);
	}
	T& getVectorValue(int *indexes)
	{
		int index=*indexes<0?(m_size+*indexes)%m_size:*indexes%m_size;
		if (m_degrees>1)
			return value.lowerVector[index]->getVectorValue(indexes+1);
		else
			return value.basicVector[index];
	}
	enum Status{Unchecked,NotInterested,Interested};
	Status checkValueAtMax(MultiVector<Status>& vec, int* indexes)
	{
		class IntValHolder
		{
		public:
			IntValHolder(int& t,int v):rel(t),val(t){rel=v;}
			~IntValHolder(){rel=val;}
		private:
			int& rel;
			int val;
		};
		T me=getVectorValue(indexes);
		auto status=vec.getVectorValue(indexes);
		if (Unchecked!=status)
			return status; 
		for (int i=0;i<m_degrees;++i)
		{
			bool kill=false;
			{
				IntValHolder holder(indexes[i],indexes[i]-1);
				if (getVectorValue(indexes)>me)
				{
					kill=true;
					checkValueAtMin(vec,indexes);
				}
			}
			if (kill)
				return vec.getVectorValue(indexes)=NotInterested;
			{
				IntValHolder holder(indexes[i],indexes[i]+1);
				if (getVectorValue(indexes)>me)
				{
					kill=true;
					checkValueAtMin(vec,indexes);
				}
			}
			if (kill)
				return vec.getVectorValue(indexes)=NotInterested;
		}
		return vec.getVectorValue(indexes)=Interested;
	}
	Status checkValueAtMin(MultiVector<Status>& vec, int* indexes)
	{
		class IntValHolder
		{
		public:
			IntValHolder(int& t,int v):rel(t),val(t){rel=v;}
			~IntValHolder(){rel=val;}
		private:
			int& rel;
			int val;
		};
		T me=getVectorValue(indexes);
		auto status=vec.getVectorValue(indexes);
		if (Unchecked!=status)
			return status; 
		for (int i=0;i<m_degrees;++i)
		{
			bool kill=false;
			{
				IntValHolder holder(indexes[i],indexes[i]-1);
				if (getVectorValue(indexes)<me)
				{
					kill=true;
					checkValueAtMin(vec,indexes);
				}
			}
			if (kill)
				return vec.getVectorValue(indexes)=NotInterested;
			{
				IntValHolder holder(indexes[i],indexes[i]+1);
				if (getVectorValue(indexes)<me)
				{
					kill=true;
					checkValueAtMin(vec,indexes);
				}
			}
			if (kill)
				return vec.getVectorValue(indexes)=NotInterested;
		}
		return vec.getVectorValue(indexes)=Interested;
	}

	int getMaxValues(std::vector<int*> &vectorIn)
	{
		int sumO=vectorIn.size();
		MultiVector<Status> o1boolVector(m_degrees,m_size);
		o1boolVector.fill(Unchecked);
		int *indexes=new int[m_degrees];
		for (int i=0;i<m_degrees;++i)
			indexes[i]=0;
		do
		if (checkValueAtMax(o1boolVector,indexes)==Interested)
			vectorIn.push_back(getIndexCopy(indexes));
		while(incressIndex(indexes));
		return vectorIn.size()-sumO;
	}
	int getMinValues(std::vector<int*> &vectorIn)
	{
		int sumO=vectorIn.size();
		MultiVector<Status> o1boolVector(m_degrees,m_size);
		o1boolVector.fill(Unchecked);
		int *indexes=new int[m_degrees];
		for (int i=0;i<m_degrees;++i)
			indexes[i]=0;
		do
		if (checkValueAtMin(o1boolVector,indexes)==Interested)
			vectorIn.push_back(getIndexCopy(indexes));
		while(incressIndex(indexes));
		return vectorIn.size()-sumO;
	}

	int getMaxValuesOld(std::vector<int*> &vectorIn)
	{
		int sumO=vectorIn.size();
		MultiVector<T>** o1recordVector=new MultiVector<T>*[m_degrees];
		MultiVector<bool>** o1boolVector=new MultiVector<bool>*[m_degrees];
		MultiVector<T>** o2recordVector=new MultiVector<T>*[m_degrees];
		for (int i=0;i<m_degrees;++i)
		{
			o1recordVector[i]=new MultiVector<T>(m_degrees,m_size);
			o1boolVector[i]=new MultiVector<bool>(m_degrees,m_size);
			o1boolVector[i]->fill(false);
			o2recordVector[i]=new MultiVector<T>(m_degrees,m_size);
		}
		int *indexes=new int[m_degrees];
		for (int i=0;i<m_degrees;++i)
			indexes[i]=0;
		for (int operatingDegree=0;operatingDegree<m_degrees;++operatingDegree)
			do
			{
				o1recordVector[operatingDegree]->getVectorValue(indexes)=get1derivative(indexes,operatingDegree);
				o2recordVector[operatingDegree]->getVectorValue(indexes)=get2derivative(indexes,operatingDegree);
			}
			while(incressIndex(indexes));
			for (int i=0;i<m_degrees;++i)
				indexes[i]=0;
			for (int operatingDegree=0;operatingDegree<m_degrees;++operatingDegree)
				do
				{
					bool tmp=o1boolVector[operatingDegree]->getVectorValue(indexes)=o1recordVector[operatingDegree]->getDiffWithNebr(indexes,operatingDegree);
					if (tmp)
						if (o2recordVector[operatingDegree]->getVectorValue(indexes)>0)
							o1boolVector[operatingDegree]->getVectorValue(indexes)=false;
				}
				while(incressIndex(indexes));

				for (int i=0;i<m_degrees;++i)
					indexes[i]=0;
				do
				{
					bool buffer=true;
					for (int operatingDegree=0;operatingDegree<m_degrees;++operatingDegree)
						if (!o1boolVector[operatingDegree]->getVectorValue(indexes))
							buffer=false;
					if (buffer)
						vectorIn.push_back(getIndexCopy(indexes));
				}
				while(incressIndex(indexes));
				return vectorIn.size()-sumO;
	}

	int getMinValuesOld(std::vector<int*> &vectorIn)
	{
		int sumO=vectorIn.size();
		MultiVector<T>** o1recordVector=new MultiVector<T>*[m_degrees];
		MultiVector<bool>** o1boolVector=new MultiVector<bool>*[m_degrees];
		MultiVector<T>** o2recordVector=new MultiVector<T>*[m_degrees];
		for (int i=0;i<m_degrees;++i)
		{
			o1recordVector[i]=new MultiVector<T>(m_degrees,m_size);
			o1boolVector[i]=new MultiVector<bool>(m_degrees,m_size);
			o1boolVector[i]->fill(false);
			o2recordVector[i]=new MultiVector<T>(m_degrees,m_size);
		}
		int *indexes=new int[m_degrees];
		for (int i=0;i<m_degrees;++i)
			indexes[i]=0;
		for (int operatingDegree=0;operatingDegree<m_degrees;++operatingDegree)
			do
			{
				o1recordVector[operatingDegree]->getVectorValue(indexes)=get1derivative(indexes,operatingDegree);
				o2recordVector[operatingDegree]->getVectorValue(indexes)=get2derivative(indexes,operatingDegree);
			}
			while(incressIndex(indexes));
			for (int i=0;i<m_degrees;++i)
				indexes[i]=0;
			for (int operatingDegree=0;operatingDegree<m_degrees;++operatingDegree)
				do
				{
					bool tmp=o1boolVector[operatingDegree]->getVectorValue(indexes)=o1recordVector[operatingDegree]->getDiffWithNebr(indexes,operatingDegree);
					if (tmp)
						if (o2recordVector[operatingDegree]->getVectorValue(indexes)<0)
							o1boolVector[operatingDegree]->getVectorValue(indexes)=false;
				}
				while(incressIndex(indexes));

				for (int i=0;i<m_degrees;++i)
					indexes[i]=0;
				do
				{
					bool buffer=true;
					for (int operatingDegree=0;operatingDegree<m_degrees;++operatingDegree)
						if (!o1boolVector[operatingDegree]->getVectorValue(indexes))
							buffer=false;
					if (buffer)
						vectorIn.push_back(getIndexCopy(indexes));
				}
				while(incressIndex(indexes));
				return vectorIn.size()-sumO;
	}
	bool getDiffWithNebr(int* index,int degreeOperating)
	{
		index[degreeOperating]+=1;
		bool v1=getVectorValue(index)>0;
		index[degreeOperating]-=2;
		bool v2=getVectorValue(index)>0;
		index[degreeOperating]+=1;
		bool v0=getVectorValue(index)>0;
		return (v1^v0)||(v2^v0);
	}
private:
	int* getIndexCopy(const int* src)
	{
		int* toReturn=new int[m_degrees];
		for (int i=0;i<m_degrees;++i)
			toReturn[i]=src[i];
		return toReturn;
	}
	double get1derivative(int* index,int degreeOperating)
	{
		index[degreeOperating]+=1;
		double v1=getVectorValue(index);
		index[degreeOperating]-=2;
		double v2=getVectorValue(index);
		index[degreeOperating]+=1;
		return v1-v2;
	}
	double get2derivative(int* index,int degreeOperating)
	{
		index[degreeOperating]+=1;
		double v1=getVectorValue(index);
		index[degreeOperating]-=2;
		double v2=getVectorValue(index);
		index[degreeOperating]+=1;
		double v0=getVectorValue(index);
		return v1+v2-v0-v0;
	}
public:
	bool incressIndex(int *indexes)
	{
		for (int i=m_degrees;i>0;--i)
		{
			indexes[i-1]++;
			if (indexes[i-1]<m_size)
				return true;
			indexes[i-1]%=m_size;
		}
		return false;
	}
};
