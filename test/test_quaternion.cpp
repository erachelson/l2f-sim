#include <iostream>
#include <vector>
#include <L2F/utils/quaternion.hpp>
#include <cassert>

void print(const L2F::quaternion &q) {
	std::cout << "(w,x,y,z) = (" << q.w << ", " << q.x << ", " << q.y << ", " << q.z << ")" << std::endl;
}

void print(const std::vector<double> &v) {
	std::cout << "(x,y,z) = (" << v.at(0) << ", " << v.at(1) << ", " << v.at(2) << ")" << std::endl;
}

void print_matrix(const std::vector<double> &m) {
	assert(m.size()==9);
	std::cout << m[0] << " " << m[1] << " " << m[2] << std::endl;
	std::cout << m[3] << " " << m[4] << " " << m[5] << std::endl;
	std::cout << m[6] << " " << m[7] << " " << m[8] << std::endl;
}

int main() {

int test=4;
/* 0 - From/to Euler
 * 1 - Rotate vector
 * 2 - Compose rotations
 * 3 - From/to matrix
 * 4 - Default constructor
 * 5 - From/to axis angle
 */
L2F::quaternion q(1., 1., 0., 0.);
L2F::quaternion p(1., 1., 0., 0.);
double pitch;
double roll;
double yaw;
std::vector<double> v;
std::vector<double> rmat;
double x, y, z, alpha;

if(test==0) {
	std::cout << "From-to Euler, should be (0.9,0.,0.)" << std::endl;
	q.fromEuler(0.9,0.,0.);
	std::cout << "Result: "; 
	q.toEuler(yaw, pitch, roll);
	std::cout << yaw << " " << pitch << " " << roll << std::endl;
	std::cout << "Angle: " << q.rotationAngle() << std::endl;
	print(q);
	std::cout << "norm = " << q.norm() << std::endl;

	std::cout << "From-to Euler, should be (0.,0.9,0.)" << std::endl;
	q.fromEuler(0.,0.9,0.);
	std::cout << "Result: "; 
	q.toEuler(yaw, pitch, roll);
	std::cout << yaw << " " << pitch << " " << roll << std::endl;
	std::cout << "Angle: " << q.rotationAngle() << std::endl;
	print(q);
	std::cout << "norm = " << q.norm() << std::endl;

	std::cout << "From-to Euler, should be (0.,0.,0.9)" << std::endl;
	q.fromEuler(0.,0.,0.9);
	std::cout << "Result: "; 
	q.toEuler(yaw, pitch, roll);
	std::cout << yaw << " " << pitch << " " << roll << std::endl;
	std::cout << "Angle: " << q.rotationAngle() << std::endl;
	print(q);
	std::cout << "norm = " << q.norm() << std::endl;

	std::cout << "From-to Euler, should be (0.,0.,-0.9)" << std::endl;
	q.fromEuler(0.,0.,-0.9);
	std::cout << "Result: "; 
	q.toEuler(yaw, pitch, roll);
	std::cout << yaw << " " << pitch << " " << roll << std::endl;
	std::cout << "Angle: " << q.rotationAngle() << std::endl;
	print(q);
	std::cout << "norm = " << q.norm() << std::endl;

	std::cout << "From-to Euler, should be (0.7,0.8,-0.9)" << std::endl;
	q.fromEuler(0.7,0.8,-0.9);
	std::cout << "Result: "; 
	q.toEuler(yaw, pitch, roll);
	std::cout << yaw << " " << pitch << " " << roll << std::endl;
	std::cout << "Angle: " << q.rotationAngle() << std::endl;
	print(q);
	std::cout << "norm = " << q.norm() << std::endl;
	
	std::cout << "From-to Euler, should be (pi/2,pi/2,pi/2)" << std::endl;
	q.fromEuler(M_PI/2.,M_PI/2.,M_PI/2.);
	std::cout << "Result: "; 
	q.toEuler(yaw, pitch, roll);
	std::cout << yaw << " " << pitch << " " << roll << std::endl;
	std::cout << "Angle: " << q.rotationAngle() << std::endl;
	print(q);
	std::cout << "norm = " << q.norm() << std::endl;
}

if(test==1) {
	v = {1., 0. ,0.};
	std::cout << "rotate of -pi/2 along yaw (Z)" << std::endl;
	print(v);
	q.fromEuler(-M_PI/2., 0., 0.);
	q.rotateVector(v);
	print(v);

	std::cout << "rotate of pi/2 along yaw (Z)" << std::endl;
	v = {1., 0. ,0.};
	print(v);
	q.fromEuler(M_PI/2., 0., 0.);
	q.rotateVector(v);
	print(v);

	std::cout << "rotate of -pi/2 along roll (X)" << std::endl;
	v = {1., 0. ,0.};
	print(v);
	q.fromEuler(0., 0., -M_PI/2.);
	q.rotateVector(v);
	print(v);

	std::cout << "rotate of pi/2 along pitch (Y)" << std::endl;
	v = {1., 0. ,0.};
	print(v);
	q.fromEuler(0., M_PI/2., 0.);
	q.rotateVector(v);
	print(v);

	std::cout << "rotate of -pi/2 along pitch (Y)" << std::endl;
	v = {1., 0. ,0.};
	print(v);
	q.fromEuler(0., -M_PI/2., 0.);
	q.rotateVector(v);
	print(v);

	std::cout << "rotate of pi/4 along yaw (Z)" << std::endl;
	v = {2., 0. ,0.};
	print(v);
	q.fromEuler(-M_PI/4., 0., 0.);
	q.rotateVector(v);
	print(v);
	
	std::cout << "rotate of pi/4 along yaw (Z), then of pi/4 along pitch." << std::endl;
	v = {1., 0. ,0.};
	print(v);
	q.fromEuler(M_PI/4., M_PI/4., 0.);
	q.rotateVector(v);
	print(v);
}

if(test==2) {
	std::cout << "compose rotations (multLeft), pi/4 along yaw then pi/4 along pitch" << std::endl;
	q.fromEuler(M_PI/4,0.,0.);
	p.fromEuler(0.,M_PI/4.,0.);
	v = {1., 0. ,0.};
	print(v);
	q.multLeft(p);
	q.rotateVector(v);
	print(v);
	std::cout << "compose rotations (multRight), pi/4 along yaw then pi/4 along pitch" << std::endl;
	q.fromEuler(M_PI/4,0.,0.);
	p.fromEuler(0.,M_PI/4.,0.);
	v = {1., 0. ,0.};
	print(v);
	p.multRight(q);
	p.rotateVector(v);
	print(v);
	std::cout << "compose rotations (multLeft), pi/4 along pitch then pi/4 along yaw" << std::endl;
	q.fromEuler(M_PI/4,0.,0.);
	p.fromEuler(0.,M_PI/4.,0.);
	v = {1., 0. ,0.};
	print(v);
	p.multLeft(q);
	p.rotateVector(v);
	print(v);
	std::cout << "compose rotations (multRight), pi/4 along pitch then pi/4 along yaw" << std::endl;
	q.fromEuler(M_PI/4,0.,0.);
	p.fromEuler(0.,M_PI/4.,0.);
	v = {1., 0. ,0.};
	print(v);
	q.multRight(p);
	q.rotateVector(v);
	print(v);
}

if(test==3) {
	std::cout << "### build quaternion for rotation of pi/2 along yaw (Z)." << std::endl;
	q.fromEuler(M_PI/2, 0., 0.);
	print(q);
	std::cout << "to matrix" << std::endl;
	q.toRotationMatrix(rmat);
	print_matrix(rmat);
	std::cout << "to quaternion" << std::endl;
	q.fromRotationMatrix(rmat);
	print(q);
	
	std::cout << "### build quaternion for rotation of pi/4 along yaw (Z), then of pi/4 along pitch." << std::endl;
	q.fromEuler(M_PI/4., M_PI/4., 0.);
	print(q);
	std::cout << "to matrix" << std::endl;
	q.toRotationMatrix(rmat);
	print_matrix(rmat);
	std::cout << "to quaternion" << std::endl;
	q.fromRotationMatrix(rmat);
	print(q);
	std::cout << "to matrix again" << std::endl;
	q.toRotationMatrix(rmat);
	print_matrix(rmat);
	
	std::cout << "### The singularity: rotations of angle pi." << std::endl;
	std::cout << "build quaternion for rotation of pi along yaw (Z)." << std::endl;
	q.fromEuler(M_PI, 0., 0.);
	print(q);
	std::cout << "to matrix" << std::endl;
	q.toRotationMatrix(rmat);
	print_matrix(rmat);
	std::cout << "to quaternion" << std::endl;
	q.fromRotationMatrix(rmat);
	print(q);
}

if(test==4) {
	L2F::quaternion qq;
	print(qq);
}

if(test==5) {
	std::cout << "Create from axis (0,0,1) and angle 0." << std::endl;
	q.fromAxisAngle(0., 0., 1., 0.);
	print(q);
	q.toAxisAngle(x, y ,z, alpha);
	std::cout << "axis = (" << x << ", " << y << ", " << z << ")" << std::endl;
	std::cout << "angle = " << 	alpha << std::endl;
	
	std::cout << "Create from axis (0,0,1) and angle pi/3=" << M_PI/3. << std::endl;
	q.fromAxisAngle(0., 0., 1., M_PI/3.);
	print(q);
	q.toAxisAngle(x, y ,z, alpha);
	std::cout << "axis = (" << x << ", " << y << ", " << z << ")" << std::endl;
	std::cout << "angle = " << 	alpha << std::endl;
	
	std::cout << "Create from axis (1,1,1) and angle -2pi/3=" << -2.*M_PI/3. << std::endl;
	q.fromAxisAngle(1., 1., 1., -2.*M_PI/3.);
	print(q);
	q.toAxisAngle(x, y ,z, alpha);
	std::cout << "axis = (" << x << ", " << y << ", " << z << ")" << std::endl;
	std::cout << "angle = " << 	alpha << std::endl;
}

}
