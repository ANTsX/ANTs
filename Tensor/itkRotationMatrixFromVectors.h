/* from
 * @article{MollerHughes99,
  author = "Tomas Mller and John F. Hughes",
  title = "Efficiently Building a Matrix to Rotate One Vector to Another",
  journal = "journal of graphics tools",
  volume = "4",
  number = "4",
  pages = "1-4",
  year = "1999",
}
http://jgt.akpeters.com/papers/MollerHughes99/code.html
*/

#include <cmath>
#include <itkVector.h>
#include <itkMatrix.h>

template <typename VectorType, typename MatrixType>
MatrixType
RotationMatrixFromVectors(VectorType from, VectorType to)
{
  double EPSILON = 0.000001;

  typename VectorType::ValueType e, h, f;
  MatrixType                     mtx;

  VectorType v = CrossProduct(from, to);

  e = from * to;
  f = (e < 0) ? -e : e;
  if (f > 1.0 - EPSILON) /* "from" and "to"-vector almost parallel */
  {
    VectorType ut, vt;     /* temporary storage vectors */
    VectorType x;          /* vector most nearly orthogonal to "from" */
    float      c1, c2, c3; /* coefficients for later use */
    int        i, j;

    x[0] = (from[0] > 0.0) ? from[0] : -from[0];
    x[1] = (from[1] > 0.0) ? from[1] : -from[1];
    x[2] = (from[2] > 0.0) ? from[2] : -from[2];

    if (x[0] < x[1])
    {
      if (x[0] < x[2])
      {
        x[0] = 1.0;
        x[1] = x[2] = 0.0;
      }
      else
      {
        x[2] = 1.0;
        x[0] = x[1] = 0.0;
      }
    }
    else
    {
      if (x[1] < x[2])
      {
        x[1] = 1.0;
        x[0] = x[2] = 0.0;
      }
      else
      {
        x[2] = 1.0;
        x[0] = x[1] = 0.0;
      }
    }

    ut[0] = x[0] - from[0];
    ut[1] = x[1] - from[1];
    ut[2] = x[2] - from[2];
    vt[0] = x[0] - to[0];
    vt[1] = x[1] - to[1];
    vt[2] = x[2] - to[2];

    c1 = 2.0 / (ut * ut);
    c2 = 2.0 / (vt * vt);
    c3 = c1 * c2 * (ut * vt);
    for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
      {
        mtx[i][j] = -c1 * ut[i] * ut[j] - c2 * vt[i] * vt[j] + c3 * vt[i] * ut[j];
      }
      mtx[i][i] += 1.0;
    }
  }
  else /* the most common case, unless "from"="to", or "from"=-"to" */
  {
    /* ...otherwise use this hand optimized version (9 mults less) */
    float hvx, hvz, hvxy, hvxz, hvyz;
    /* h = (1.0 - e)/DOT(v, v); old code */
    h = 1.0 / (1.0 + e); /* optimization by Gottfried Chen */
    hvx = h * v[0];
    hvz = h * v[2];
    hvxy = hvx * v[1];
    hvxz = hvx * v[2];
    hvyz = hvz * v[1];
    mtx[0][0] = e + hvx * v[0];
    mtx[0][1] = hvxy - v[2];
    mtx[0][2] = hvxz + v[1];

    mtx[1][0] = hvxy + v[2];
    mtx[1][1] = e + h * v[1] * v[1];
    mtx[1][2] = hvyz - v[0];

    mtx[2][0] = hvxz - v[1];
    mtx[2][1] = hvyz + v[0];
    mtx[2][2] = e + hvz * v[2];
  }

  return mtx;
}

// /*
//  * A function for creating a rotation matrix that rotates a vector called
//  * "from" into another vector called "to".
//  * Input : from[3], to[3] which both must be *normalized* non-zero vectors
//  * Output: mtx[3][3] -- a 3x3 matrix in colum-major form
//  * Authors: Tomas Möller, John Hughes
//  *          "Efficiently Building a Matrix to Rotate One Vector to Another"
//  *          Journal of Graphics Tools, 4(4):1-4, 1999
//  */
// void fromToRotation(float from[3], float to[3], float mtx[3][3]) {
//   float v[3];
//   float e, h, f;

//   CROSS(v, from, to);
//   e = DOT(from, to);
//   f = (e < 0)? -e:e;
//   if (f > 1.0 - EPSILON)     /* "from" and "to"-vector almost parallel */
//   {
//     float u[3], v[3]; /* temporary storage vectors */
//     float x[3];       /* vector most nearly orthogonal to "from" */
//     float c1, c2, c3; /* coefficients for later use */
//     int i, j;

//     x[0] = (from[0] > 0.0)? from[0] : -from[0];
//     x[1] = (from[1] > 0.0)? from[1] : -from[1];
//     x[2] = (from[2] > 0.0)? from[2] : -from[2];

//     if (x[0] < x[1])
//     {
//       if (x[0] < x[2])
//       {
//         x[0] = 1.0; x[1] = x[2] = 0.0;
//       }
//       else
//       {
//         x[2] = 1.0; x[0] = x[1] = 0.0;
//       }
//     }
//     else
//     {
//       if (x[1] < x[2])
//       {
//         x[1] = 1.0; x[0] = x[2] = 0.0;
//       }
//       else
//       {
//         x[2] = 1.0; x[0] = x[1] = 0.0;
//       }
//     }

//     u[0] = x[0] - from[0]; u[1] = x[1] - from[1]; u[2] = x[2] - from[2];
//     v[0] = x[0] - to[0];   v[1] = x[1] - to[1];   v[2] = x[2] - to[2];

//     c1 = 2.0 / DOT(u, u);
//     c2 = 2.0 / DOT(v, v);
//     c3 = c1 * c2  * DOT(u, v);

//     for (i = 0; i < 3; i++) {
//       for (j = 0; j < 3; j++) {
//         mtx[i][j] =  - c1 * u[i] * u[j]
//                      - c2 * v[i] * v[j]
//                      + c3 * v[i] * u[j];
//       }
//       mtx[i][i] += 1.0;
//     }
//   }
//   else  /* the most common case, unless "from"="to", or "from"=-"to" */
//   {
// #if 0
//     /* unoptimized version - a good compiler will optimize this. */
//     /* h = (1.0 - e)/DOT(v, v); old code */
//     h = 1.0/(1.0 + e);      /* optimization by Gottfried Chen */
//     mtx[0][0] = e + h * v[0] * v[0];
//     mtx[0][1] = h * v[0] * v[1] - v[2];
//     mtx[0][2] = h * v[0] * v[2] + v[1];

//     mtx[1][0] = h * v[0] * v[1] + v[2];
//     mtx[1][1] = e + h * v[1] * v[1];
//     mtx[1][2] = h * v[1] * v[2] - v[0];

//     mtx[2][0] = h * v[0] * v[2] - v[1];
//     mtx[2][1] = h * v[1] * v[2] + v[0];
//     mtx[2][2] = e + h * v[2] * v[2];
// #else
//     /* ...otherwise use this hand optimized version (9 mults less) */
//     float hvx, hvz, hvxy, hvxz, hvyz;
//     /* h = (1.0 - e)/DOT(v, v); old code */
//     h = 1.0/(1.0 + e);      /* optimization by Gottfried Chen */
//     hvx = h * v[0];
//     hvz = h * v[2];
//     hvxy = hvx * v[1];
//     hvxz = hvx * v[2];
//     hvyz = hvz * v[1];
//     mtx[0][0] = e + hvx * v[0];
//     mtx[0][1] = hvxy - v[2];
//     mtx[0][2] = hvxz + v[1];

//     mtx[1][0] = hvxy + v[2];
//     mtx[1][1] = e + h * v[1] * v[1];
//     mtx[1][2] = hvyz - v[0];

//     mtx[2][0] = hvxz - v[1];
//     mtx[2][1] = hvyz + v[0];
//     mtx[2][2] = e + hvz * v[2];
// #endif
//   }
// }
