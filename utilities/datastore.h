/**
* Datastructure for storing our cell data in 3d vectors (for positions, velocities etc) and in 4d (particle data in flex)
*
**/
#ifndef DATASTORE_H_
#define DATASTORE_H_

#include<cmath>
#include<string>


	/** 3d integer vector */
	struct intVec3d {
		int x, y, z;
	};

	/** 3d double vector */
	struct doubleVec3d {
		double x, y, z;
	};

	/** 3d float vector */
	struct floatVec3d {
		float x, y, z;
	};

	/** 4d float vector (useful for FleX particle data (x,y,z,1/m))*/
	struct floatVec4d {
		float x, y, z, m;
	};


	class myVec3d {
	public:
		myVec3d();
		myVec3d(double x, double y, double z);
		myVec3d(const myVec3d& v);
		~myVec3d();
		myVec3d& operator=(myVec3d rhs) {
			pos.x = rhs.pos.x;
			pos.y = rhs.pos.y;
			pos.z = rhs.pos.z;
			return *this;
		}
		;
		void normalize();
		double length();
		double lengthSquared();
		void add(const myVec3d& v);
		void sub(const myVec3d& v);
		void negate(const myVec3d& v);
		void scale(double s);
		void scale(double s, const myVec3d& v);
		void cross(const myVec3d& v1, const myVec3d& v2);
		doubleVec3d pos;
	};

	template<typename T>
	class CellarisVector {
	private:
		T * elem;
		unsigned long size;
		unsigned long count;
		unsigned long capacity;
	public:
		unsigned long getSize() {
			return size;
		}
		unsigned long getCount() {
			return count;
		}
		CellarisVector() {
			elem = new T[10];
			count = size = 0;
			capacity = 10;
		}
		CellarisVector(unsigned long s) {
			elem = new T[s];
			count = size = 0;
			capacity = s;
		}
		~CellarisVector() {
			delete[] elem;
		}
		void setCapacity(unsigned long s) {
			elem = new T[s];
			count = size = 0;
			capacity = s;
		}
		T& operator[](int i) {
			return elem[i];
		}
		void add(T& newElem) {
			elem[size++] = newElem;
			count++;
			if (size / capacity > 0.8) {
				capacity *= 2;
				T* temp = new T[capacity];
				int j = 0;
				for (int i = 0; i != size; ++i) {
					if (elem[i] != NULL)
						temp[j++] = elem[i];
				}
				delete[] elem;
				elem = temp;
				count = size = j;
			}
		}
		;
		void remove(T& newElem) {
			for (int i = 0; i != size; ++i)
				if (newElem == elem[i]) {
					elem[i] = NULL;
					count--;
					return;
				}
		}
		;
	};

#endif /* DATASTOR_H_ */