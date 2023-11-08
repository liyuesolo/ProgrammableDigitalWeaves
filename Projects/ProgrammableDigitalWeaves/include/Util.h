#ifndef UTIL_H
#define UTIL_H

#include "VecMatDef.h"

//http://paulbourke.net/geometry/circlesphere/tvoght.c
template<class T, int dim>
bool circleCircleIntersection(const Vector<T, dim>& x0, T r0,
                               const Vector<T, dim>& x1, T r1,
                               Vector<T, dim>& ixn0,
                               Vector<T, dim>& ixn1)
{
    ixn0 = Vector<T, dim>::Zero();
    ixn1 = Vector<T, dim>::Zero();
    T a, dx, dy, d, h, rx, ry;
    T x2, y2;

    /* dx and dy are the vertical and horizontal distances between
    * the circle centers.
    */
    dx = x1[0] - x0[0];
    dy = x1[1] - x0[1];

    /* Determine the straight-line distance between the centers. */
    //d = sqrt((dy*dy) + (dx*dx));
    d = std::hypot(dx,dy); // Suggested by Keith Briggs

    /* Check for solvability. */
    if (d > (r0 + r1))
    {
    /* no solution. circles do not intersect. */
    return false;
    }
    if (d < fabs(r0 - r1))
    {
    /* no solution. one circle is contained in the other */
    return false;
    }

    /* 'point 2' is the point where the line through the circle
    * intersection points crosses the line between the circle
    * centers.  
    */

    /* Determine the distance from point 0 to point 2. */
    a = ((r0*r0) - (r1*r1) + (d*d)) / (2.0 * d) ;

    /* Determine the coordinates of point 2. */
    x2 = x0[0] + (dx * a/d);
    y2 = x0[1] + (dy * a/d);

    /* Determine the distance from point 2 to either of the
    * intersection points.
    */
    h = std::sqrt((r0*r0) - (a*a));

    /* Now determine the offsets of the intersection points from
    * point 2.
    */
    rx = -dy * (h/d);
    ry = dx * (h/d);

    /* Determine the absolute intersection points. */
    ixn0[0] = x2 + rx;
    ixn1[0] = x2 - rx;
    ixn0[1] = y2 + ry;
    ixn1[1] = y2 - ry;

    return true;
}

template<class T>
inline T cross2D(const Vector<T, 2> x1, const Vector<T, 2> x2)
{
	 return x1[0] * x2[1] - x1[1]*x2[0];
}

//https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
template<class T>
bool lineSegementsIntersect2D(const Vector<T, 2>& p, const Vector<T, 2>& p2, const Vector<T, 2>& q, const Vector<T, 2>& q2, 
    Vector<T, 2>& intersection, bool considerCollinearOverlapAsIntersect = true)
{
  using TV = Vector<T, 2>;

  intersection = TV::Zero();

  TV r = p2 - p;
  TV s = q2 - q;
  TV qmp = q - p;

  T rxs = cross2D(r, s);
  T qpxr = cross2D(qmp, r);

  if (std::abs(rxs) < 1e-6 && std::abs(qpxr) < 1e-6)
  {
      // 1. If either  0 <= (q - p) * r <= r * r or 0 <= (p - q) * s <= * s
      // then the two lines are overlapping,
      if (considerCollinearOverlapAsIntersect)
          if ((0 <= (q - p).dot(r) && (q - p).dot(r) <= r.dot(r)) || (0 <= (p - q).dot(s) && (p - q).dot(s) <= s.dot(s)))
              return true;

      // 2. If neither 0 <= (q - p) * r = r * r nor 0 <= (p - q) * s <= s * s
      // then the two lines are collinear but disjoint.
      // No need to implement this expression, as it follows from the expression above.
      return false;
  }

  // 3. If r x s = 0 and (q - p) x r != 0, then the two lines are parallel and non-intersecting.
  if (std::abs(rxs) < 1e-6 && std::abs(qpxr) > 1e-6)
      return false;

  // t = (q - p) x s / (r x s)
  T t = cross2D(qmp, s)/rxs;

  // u = (q - p) x r / (r x s)

  T u = cross2D(qmp, r) / rxs;

  // 4. If r x s != 0 and 0 <= t <= 1 and 0 <= u <= 1
  // the two line segments meet at the point p + t r = q + u s.
  if (std::abs(rxs) > 1e-6 && (0 <= t && t <= 1) && (0 <= u && u <= 1))
  {
      // We can calculate the intersection point using either t or u.
      intersection = p + t*r;
      // An intersection was found.
      return true;
  }

  // 5. Otherwise, the two line segments are not parallel but do not intersect.
  return false;
}


template<class T>
inline double Det(const Vector<T, 2> x1, const Vector<T, 2> x2)
{
	return x1[0] * x2[1] - x1[1]*x2[0];
}

///Calculate intersection of two lines.
///\return true if found, false if not found or error
template<class T>
bool lineLineIntersect(const Vector<T, 2>& x1, //Line 1 start
	const Vector<T, 2>& x2, //Line 1 end
	const Vector<T, 2>& x3, //Line 2 start
	const Vector<T, 2>& x4, //Line 2 end
	Vector<T, 2>& intersection) //Output 
{
	//http://mathworld.wolfram.com/Line-LineIntersection.html

	T detL1 = Det(x1, x2);
	T detL2 = Det(x3, x4);
	T x1mx2 = x1[0] - x2[0];
	T x3mx4 = x3[0] - x4[0];
	T y1my2 = x1[1] - x2[1];
	T y3my4 = x3[1] - x4[1];

	T xnom = Det(Vector<T, 2>(detL1, x1mx2), Vector<T, 2>(detL2, x3mx4));
	T ynom = Det(Vector<T, 2>(detL1, y1my2), Vector<T, 2>(detL2, y3my4));
	T denom = Det(Vector<T, 2>(x1mx2, y1my2), Vector<T, 2>(x3mx4, y3my4));
	if(denom == 0.0)//Lines don't seem to cross
	{
		intersection = Vector<T, 2>(NAN, NAN);
		return false;
	}

	intersection[0] = xnom / denom;	
	intersection[1] = ynom / denom;
	if(!std::isfinite(intersection[0]) || !std::isfinite(intersection[1])) //Probably a numerical issue
		return false;

	return true; //All OK
}

//https://www.codeproject.com/Tips/84226/Is-a-Point-inside-a-Polygon
template<class T>
bool insidePolygon(const std::vector<Vector<T, 2>>& parallogram, const Vector<T, 2>& test)
{
    int nvert = parallogram.size();
    int i, j;
    bool c = false;
    for (i = 0, j = nvert-1; i < nvert; j = i++) {
        if ( ((parallogram[i][1]>test[1]) != (parallogram[j][1]>test[1])) &&
        (test[0] < (parallogram[j][0]-parallogram[i][0]) * (test[1]-parallogram[i][1]) / (parallogram[j][1]-parallogram[i][1]) + parallogram[i][0]) )
        c = !c;
    }
    return c;
}

template<class VectorType>
inline void appendToEnd(std::vector<VectorType>& a, const std::vector<VectorType>& b)
{
    a.insert(a.end(), b.begin(), b.end());
}

template<class T>
inline T signedAngle(const Vector<T, 3>& u, const Vector<T, 3>& v, const Vector<T, 3>& n) {
  Vector<T, 3> w = u.cross(v);
  T angle = std::atan2(w.norm(), u.dot(v));
  if (n.dot(w) < 0) return -angle;
  return angle;
}

template<class T, int dim>
inline bool colinear(Vector<T, dim> a, Vector<T, dim> b)
{
    if((a-b).norm()<1e-2)
		return true;
	if((a-b).norm()>1.99)
		return true;
	return false;
}

template <class T>
inline void rotateAxisAngle(Vector<T, 3>& v,
                            const Vector<T, 3>& z,
                            const T theta) 
{
  
  if (theta == 0) return;

  T c = cos(theta);
  T s = sin(theta);

  v = c * v + s * z.cross(v) + z.dot(v) * (1.0 - c) * z;
}

template<class T>
Matrix<T, 3, 3> rotationMatrixFromEulerAngle(T angle_z, T angle_y, T angle_x)
{
  Matrix<T, 3, 3> R, yaw, pitch, roll;
  yaw.setZero(); pitch.setZero(); roll.setZero();
  yaw(0, 0) = cos(angle_z);	yaw(0, 1) = -sin(angle_z);
  yaw(1, 0) = sin(angle_z);	yaw(1, 1) = cos(angle_z);
  yaw(2, 2) = 1.0;
  //y rotation
  pitch(0, 0) = cos(angle_y); pitch(0, 2) = sin(angle_y);
  pitch(1, 1) = 1.0;
  pitch(2, 0) = -sin(angle_y); pitch(2, 2) = cos(angle_y);
  //x rotation
  roll(0, 0) = 1.0;
  roll(1, 1) = cos(angle_x); roll(1, 2) = -sin(angle_x);
  roll(2, 1) = sin(angle_x); roll(2, 2) = cos(angle_x);
  R = yaw * pitch * roll;
  return R;
}

template<class T>
Vector<T, 3> parallelTransport(const Vector<T, 3>& u, const Vector<T, 3>& t0, const Vector<T, 3>& t1) 
{
  
    Vector<T, 3> b = t0.cross(t1);


    if(b.norm() < std::numeric_limits<T>::epsilon())
        return u;

    b.normalize();

    Vector<T, 3> n0 = t0.cross(b).normalized();
    Vector<T, 3> n1 = t1.cross(b).normalized();

    return u.dot(t0.normalized()) * t1.normalized() + u.dot(n0) * n1 +
            u.dot(b) * b;
}

template<class T>
Vector<T, 3> parallelTransportOrthonormalVector(const Vector<T, 3>& u, const Vector<T, 3>& t0, const Vector<T, 3>& t1) {
  
    Vector<T, 3> b = t0.cross(t1);


    if(b.norm() < std::numeric_limits<T>::epsilon())
        return u;

    b.normalize();

    Vector<T, 3> n0 = t0.cross(b);
    Vector<T, 3> n1 = t1.cross(b);

    return u.dot(n0) * n1 + u.dot(b) * b;
}

#endif