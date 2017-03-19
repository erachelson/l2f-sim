#ifndef L2FSIM_QUATERNION_HPP_
#define L2FSIM_QUATERNION_HPP_

#include <cmath>
#include <vector>
#include <algorithm>

namespace L2Fsim {

/**
 * Quaternion class for the L2F library
 * Convention (for Euler angles - the rest is independent of this convention):
 * - yaw is a rotation around the z-axis
 * - pitch is a rotation around the y-axis
 * - roll is a rotation around the x-axis
 * Reference: Euler Angles, Quaternions, and Transformation Matrices.
 * NASA-TM-74839, shuttle program (1977).
 */
class quaternion {
protected:
	double w;
	double x;
	double y;
	double z;

public:
	/** Constructor (default is the identity rotation). */
	quaternion(double w_=1., double x_=0., double y_=0., double z_=0.);
	/** Copy constructor */
	quaternion(const quaternion &q2);
	/** Destructor */
	virtual ~quaternion();
	/** Initializes from Euler angles */
	void fromEuler(double yaw, double pitch, double roll);
	/** Returns corresponding Euler angles */
	void toEuler(double &yaw, double &pitch, double &roll) const;
	/** Initializes from a pair axis-angle. Axis should be a non-zero vector (does not need to be normalized). */
	void fromAxisAngle(double xx, double yy, double zz, double alpha);
	/** Returns the corresponding (not normalized) rotation axis and angle. If angle is zero, then the returned axis is a zero vector. */
	void toAxisAngle(double &xx, double &yy, double &zz, double &alpha) const;
	/** Initializes from a rotation matrix. The argument is supposed to be a 9-elements vector containing the consecutive rows of the matrix. This function supposes that the provided matrix is indeed a rotation matrix and does not perform any verification on it. */
	void fromRotationMatrix(const std::vector<double>& m);
	/** Returns the corresponding rotation matrix. The return value is a 9-elements vector containing the consecutive rows of the matrix. */
	void toRotationMatrix(std::vector<double> &m) const;
	/** Rotates vector v using the quaternion */
	void rotateVector(std::vector<double> &v) const;
	/** Changes q into q*q2 (so the corresponding rotation becomes "r2 then r"). */
	void multRight(quaternion &q2);
	/** Changes q into q2*q (so the corresponding rotation becomes "r then r2"). */
	void multLeft(quaternion &q2);
	/** Returns the rotation angle */
	double rotationAngle() const;
	/** Returns the (not normalized) rotation axis */
	void rotationAxis(std::vector<double> &v) const;
	/** Returns the (not normalized) rotation axis */
	void rotationAxis(double &xx, double &yy, double &zz) const;
	/** Normalizes the quaternion */
	void normalize();
	/** Reverses the rotation angle */
	void invert();
	/** The quaternion's norm */
	double norm() const;
};

quaternion::quaternion(double w_, double x_, double y_, double z_)
: w(w_), x(x_), y(y_), z(z_) { this->normalize(); }

quaternion::quaternion(const quaternion &q2)
: w(q2.w), x(q2.x), y(q2.y), z(q2.z) { }

quaternion::~quaternion() { }

void quaternion::fromEuler(double yaw, double pitch, double roll) {
	double c1 = std::cos(yaw / 2.);
	double s1 = std::sin(yaw / 2.);
	double c2 = std::cos(pitch / 2.);
	double s2 = std::sin(pitch / 2.);
	double c3 = std::cos(roll / 2.);
	double s3 = std::sin(roll / 2.);
	w = c1*c2*c3 + s1*s2*s3;
	x = c1*c2*s3 - s1*s2*c3;
	y = c1*s2*c3 + s1*c2*s3;
	z = s1*c2*c3 - c1*s2*s3;
}

void quaternion::toEuler(double &yaw, double &pitch, double &roll) const {
	double sqw = w*w;
	double sqx = x*x;
	double sqy = y*y;
	double sqz = z*z;
	yaw   = std::atan2(2.*(x*y + z*w), sqw+sqx-sqy-sqz);
	roll  = std::atan2(2.0 * (y*z + x*w),(-sqx - sqy + sqz + sqw));
	pitch = std::asin(-2.0 * (x*z - y*w)/(sqx + sqy + sqz + sqw));
}

void quaternion::fromAxisAngle(double xx, double yy, double zz, double alpha) {
	double norm = std::sqrt(xx*xx + yy*yy + zz*zz);
	w = std::cos(alpha/2.);
	double s = std::sin(alpha/2.);
	x = xx*s/norm;
	y = yy*s/norm;
	z = zz*s/norm;
}

void quaternion::toAxisAngle(double &xx, double &yy, double &zz, double &alpha) const {
	alpha = 2.*std::acos(w);
	xx = x;
	yy = y;
	zz = z;
}

void quaternion::fromRotationMatrix(const std::vector<double>& m) {
	w = std::sqrt( std::max( 0., 1. + m[0] + m[4] + m[8] ) ) / 2.;
	x = std::sqrt( std::max( 0., 1. + m[0] - m[4] - m[8] ) ) / 2.;
	y = std::sqrt( std::max( 0., 1. - m[0] + m[4] - m[8] ) ) / 2.;
	z = std::sqrt( std::max( 0., 1. - m[0] - m[4] + m[8] ) ) / 2.;
	x = std::copysign(x, m[7] - m[5]);
	y = std::copysign(y, m[2] - m[6]);
	z = std::copysign(z, m[3] - m[1]);
	this->normalize();
}

void quaternion::toRotationMatrix(std::vector<double> &m) const {
	std::vector<double>(9).swap(m);
	double wx = w*x;
	double wy = w*y;
	double wz = w*z;
	double xx = x*x;
	double xy = x*y;
	double xz = x*z;
	double yy = y*y;
	double yz = y*z;
	double zz = z*z;
	/* 1st row */
	m[0] = 1.-2.*(yy+zz);
	m[1] =    2.*(xy-wz);
	m[2] =    2.*(wy+xz);
	/* 2nd row */
	m[3] =    2.*(wz+xy);
	m[4] = 1.-2.*(xx+zz);
	m[5] =    2.*(yz-wx);
	/* 3rd row */
	m[6] =    2.*(xz-wy);
	m[7] =    2.*(wx+yz);
	m[8] = 1.-2.*(xx+yy);

}

void quaternion::rotateVector(std::vector<double> &v) const {
	double wx = w*x;
	double wy = w*y;
	double wz = w*z;
	double xx = x*x;
	double xy = x*y;
	double xz = x*z;
	double yy = y*y;
	double yz = y*z;
	double zz = z*z;
	double v1 = 2.*( (-yy - zz)*v[0] + ( xy - wz)*v[1] + ( wy + xz)*v[2] ) + v[0];
	double v2 = 2.*( ( wz + xy)*v[0] + (-xx - zz)*v[1] + ( yz - wx)*v[2] ) + v[1];
	double v3 = 2.*( ( xz - wy)*v[0] + ( wx + yz)*v[1] + (-xx - yy)*v[2] ) + v[2];
	std::vector<double>({v1,v2,v3}).swap(v);
}

void quaternion::multRight(quaternion &q2) {
	double ww = w*q2.w - x*q2.x - y*q2.y - z*q2.z;
	double xx = w*q2.x + x*q2.w + y*q2.z - z*q2.y;
	double yy = w*q2.y + y*q2.w + z*q2.x - x*q2.z;
	double zz = w*q2.z + z*q2.w + x*q2.y - y*q2.x;
	w = ww;
	x = xx;
	y = yy;
	z = zz;
}

void quaternion::multLeft(quaternion &q2) {
	double ww = q2.w*w - q2.x*x - q2.y*y - q2.z*z;
	double xx = q2.w*x + q2.x*w + q2.y*z - q2.z*y;
	double yy = q2.w*y + q2.y*w + q2.z*x - q2.x*z;
	double zz = q2.w*z + q2.z*w + q2.x*y - q2.y*x;
	w = ww;
	x = xx;
	y = yy;
	z = zz;

}

double quaternion::rotationAngle() const { return 2.*std::acos(w); }

void quaternion::rotationAxis(std::vector<double> &v) const {
	std::vector<double>({x,y,z}).swap(v);
}

void quaternion::rotationAxis(double &xx, double &yy, double &zz) const {
	xx = x; yy = y; zz = z;
}

void quaternion::normalize() {
	double magnitude = sqrt(w*w + x*x + y*y + z*z);
	w /= magnitude;
	x /= magnitude;
	y /= magnitude;
	z /= magnitude;
}

void quaternion::invert() {
	x = -x;
	y = -y;
	z = -z;
}

double quaternion::norm() const {
	return std::sqrt(w*w+x*x+y*y+z*z);
}

}

#endif
