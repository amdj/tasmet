// Adding mode required functions over here...


inline void printVec(vd vec) {
	us i;
	cout << "[" ;
	for (i=0;i<vec.size()-1;i++)  {
	cout << vec[i] << " ";
	}
	cout << vec[vec.size()-1] << "]" << endl;
}

// We add what we need...
//inline vd operator+(const vd &v1,const vd &v2) {
//	assert(v1.size()==v2.size());
//	vd result=vd(v1.size());
//	unsigned i;
//	for (i=0;i<v1.size() ;i++)
//		result[i]=(v1)[i]+(v2)[i];
//	return result;
//}
//inline vd operator-(const vd &v1,const vd &v2) {
//	assert(v1.size()==v2.size());
//	vd result=vd(v1.size());
//	unsigned i;
//	for (i=0;i<v1.size() ;i++)
//		result[i]=(v1)[i]-(v2)[i]; // Very subtle minus here..
//	return result;
//}
//inline vd operator/(const vd &v1,double v2) {
//	vd result=vd(v1.size());
//	unsigned i;
//	for (i=0;i<v1.size() ;i++)
//		result[i]=(v1)[i]/v2;
//	return result;
//}
//inline vd operator*(const vd &v1,const vd &v2) {
//	assert(v1.size()==v2.size());
//	vd result=vd(v1.size());
//	unsigned i;
//	for (i=0;i<v1.size() ;i++)
//		result[i]=v1[i]*v2[i]; // Very subtle minus here..
//	return result;
//}
//inline vd operator*(double v1,const vd &v2) {
//	vd result=vd(v2.size());
//	us i;
//	for (i=0;i<v2.size() ;i++)
//		result[i]=v2[i]*v1;
//	return result;
//}
//inline vd operator*(const vd &v1,double v2) {
//	//assert(v1.size()==v2.size());
//	vd result=vd(v1.size());
//	unsigned i;
//	for (i=0;i<v1.size() ;i++)
//		result[i]=v1[i]*v2; // Very subtle minus here..
//	return result;
//}
//inline vd operator^(const vd &v1,double v2) {
//	vd result=vd(v1.size());
//	unsigned i;
//	for (i=0;i<v1.size() ;i++)
//		result[i]=pow(v1[i],v2); // Very subtle minus here..
//	return result;
//}

