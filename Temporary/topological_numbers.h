#ifndef TOPO_INCLUDED
#define TOPO_INCLUDED
/* This file provides some code to compute the topological numbers
   introduced by Giles Bertrand and Gregoire Malandin.

   This code is the result of my understanding of their paper, and is
   not optimized. In fact, it is possible to generate look-up tables
   to greatly speed-up computations...

   Given a binary 3*3*3 neighborhood ( the structure
   TOPOLOGICAL_NEIGHBORHOOD or NBH ) and a pair of consistent
   connectivities ( 1=(6+,18), 2=(18,6+), 3=(6,26), 4=(26,6) ),
   the function checkTn computes the topological number associated with
   the structure NBH and a specific connectivity. The fucntion
   checkSimple checks if the point is simple or not (see definitions below).

   The topological neighborhood NBH has to be initialized as a binary
   object 0-1, where 1 represents the Foreground Object and 0 the
   Background Object. The first connectivity in the pair of digital
   connectivities ( 1=(6+,18), 2=(18,6+), 3=(6,26), 4=(26,6) ) refers to
   the Foreground Object, the second one to the Background Object
*/

/*
//////////////////////////////////////////////////////////////////////
//     COMPATIBLE CONNECTIVITIES / TOPOLOGICAL NUMBERS
//
//     TOPOLOGICAL CONVENTION
//     0:            No topological constraint
//     1:            (6+,18)
//     2:            (18,6+)
//     3:            (6,26)
//     4:            (26,6)
//     default:      (6+,18)
////////////////////////////////////////////////////////////////////
*/

typedef unsigned char TOPOLOGICAL_NEIGHBORHOOD[3][3][3];

#define NBH TOPOLOGICAL_NEIGHBORHOOD
#define TMP 100
#define MAXIMUM_NUMBER_OF_COMPONENTS 10 // maximum number of components in a NBH
#define MAX_COMP MAXIMUM_NUMBER_OF_COMPONENTS

/*The connectivity number is equal to the maximum allowed distance
  (||.||1 norm) used to generate the topological neighborhood*/
int
connectivityNumber(int connectivity)
{
  int con;

  switch (connectivity)
  {
    case 1: /*for the connectivity 1 [=(6+,18)], the topological
      neighborhood of a point x is defined by {x /
      ||a-x||1<=1 } -> the 6 neighborhood! */
      con = 1;
      break;
    case 2: /*the 18-neighborhood*/
      con = 2;
      break;
    case 3: /*the 6-heighborhood*/
      con = 1;
      break;
    case 4: /*the 26-neighborhood*/
      con = 3;
      break;
    default:
      con = 1;
      break;
  }
  return con;
}

// Digital topology requires a pair of compatible connectivities
int
associatedConnectivity(int connectivity)
{
  switch (connectivity)
  {
    case 1:
      return 2;
      break;
    case 2:
      return 1;
      break;
    case 3:
      return 4;
      break;
    case 4:
      return 3;
      break;
    default:
      return 2;
      break;
  }
  return 0;
}

NBH *
reverseNBH(NBH * nbh_src, NBH * nbh_dst)
{
  NBH * nbh;
  int   i, j, k;

  if (nbh_dst)
  {
    nbh = nbh_dst;
  }
  else
  {
    nbh = (NBH *)malloc(sizeof(NBH));
  }
  for (k = 0; k < 3; k++)
  {
    for (j = 0; j < 3; j++)
    {
      for (i = 0; i < 3; i++)
      {
        if ((*nbh_src)[i][j][k])
        {
          (*nbh)[i][j][k] = 0;
        }
        else
        {
          (*nbh)[i][j][k] = 1;
        }
      }
    }
  }

  return nbh;
}

NBH *
N_6_1(NBH * nbh_src, NBH * nbh_dst)
{
  int   i, j, k;
  NBH * nbh;

  if (!nbh_dst)
  {
    nbh = (NBH *)calloc(1, sizeof(NBH));
    if ((*nbh_src)[0][1][1])
    {
      (*nbh)[0][1][1] = 1;
    }
    if ((*nbh_src)[2][1][1])
    {
      (*nbh)[2][1][1] = 1;
    }
    if ((*nbh_src)[1][0][1])
    {
      (*nbh)[1][0][1] = 1;
    }
    if ((*nbh_src)[1][2][1])
    {
      (*nbh)[1][2][1] = 1;
    }
    if ((*nbh_src)[1][1][0])
    {
      (*nbh)[1][1][0] = 1;
    }
    if ((*nbh_src)[1][1][2])
    {
      (*nbh)[1][1][2] = 1;
    }
    return nbh;
  }
  else
  {
    nbh = nbh_dst;
  }

  if ((*nbh_src)[0][1][1])
  {
    (*nbh)[0][1][1] = TMP;
  }
  if ((*nbh_src)[2][1][1])
  {
    (*nbh)[2][1][1] = TMP;
  }
  if ((*nbh_src)[1][0][1])
  {
    (*nbh)[1][0][1] = TMP;
  }
  if ((*nbh_src)[1][2][1])
  {
    (*nbh)[1][2][1] = TMP;
  }
  if ((*nbh_src)[1][1][0])
  {
    (*nbh)[1][1][0] = TMP;
  }
  if ((*nbh_src)[1][1][2])
  {
    (*nbh)[1][1][2] = TMP;
  }
  for (k = 0; k < 3; k++)
  {
    for (j = 0; j < 3; j++)
    {
      for (i = 0; i < 3; i++)
      {
        if ((*nbh)[i][j][k] == TMP)
        {
          (*nbh)[i][j][k] = 1;
        }
        else
        {
          (*nbh)[i][j][k] = 0;
        }
      }
    }
  }

  return nbh;
}

NBH *
N_6_2(NBH * nbh_src, NBH * nbh_dst)
{
  int   i, j, k;
  NBH * nbh;

  nbh = N_6_1(nbh_src, nullptr);
  for (i = 0; i < 3; i = i + 2)
  {
    if ((*nbh)[i][1][1])
    {
      if ((*nbh_src)[i][0][1])
      {
        (*nbh)[i][0][1] = 1;
      }
      if ((*nbh_src)[i][2][1])
      {
        (*nbh)[i][2][1] = 1;
      }
      if ((*nbh_src)[i][1][0])
      {
        (*nbh)[i][1][0] = 1;
      }
      if ((*nbh_src)[i][1][2])
      {
        (*nbh)[i][1][2] = 1;
      }
    }
  }
  for (j = 0; j < 3; j = j + 2)
  {
    if ((*nbh)[1][j][1])
    {
      if ((*nbh_src)[0][j][1])
      {
        (*nbh)[0][j][1] = 1;
      }
      if ((*nbh_src)[2][j][1])
      {
        (*nbh)[2][j][1] = 1;
      }
      if ((*nbh_src)[1][j][0])
      {
        (*nbh)[1][j][0] = 1;
      }
      if ((*nbh_src)[1][j][2])
      {
        (*nbh)[1][j][2] = 1;
      }
    }
  }
  for (k = 0; k < 3; k = k + 2)
  {
    if ((*nbh)[1][1][k])
    {
      if ((*nbh_src)[0][1][k])
      {
        (*nbh)[0][1][k] = 1;
      }
      if ((*nbh_src)[2][1][k])
      {
        (*nbh)[2][1][k] = 1;
      }
      if ((*nbh_src)[1][0][k])
      {
        (*nbh)[1][0][k] = 1;
      }
      if ((*nbh_src)[1][2][k])
      {
        (*nbh)[1][2][k] = 1;
      }
    }
  }

  if (nbh_dst)
  {
    for (k = 0; k < 3; k++)
    {
      for (j = 0; j < 3; j++)
      {
        for (i = 0; i < 3; i++)
        {
          (*nbh_dst)[i][j][k] = (*nbh)[i][j][k];
        }
      }
    }
    free(nbh);
    nbh = nbh_dst;
  }

  return nbh;
}

NBH *
N_6_3(NBH * nbh_src, NBH * nbh_dst)
{
  int   i, j, k;
  NBH * nbh;

  nbh = N_6_2(nbh_src, nullptr);

  i = 0;
  j = 0;
  k = 0;
  if ((*nbh_src)[i][j][k])
  {
    if ((*nbh)[1][j][k] || (*nbh)[i][1][k] || (*nbh)[i][j][1])
    {
      (*nbh)[i][j][k] = 1;
    }
  }
  i = 0;
  j = 0;
  k = 2;
  if ((*nbh_src)[i][j][k])
  {
    if ((*nbh)[1][j][k] || (*nbh)[i][1][k] || (*nbh)[i][j][1])
    {
      (*nbh)[i][j][k] = 1;
    }
  }
  i = 0;
  j = 2;
  k = 0;
  if ((*nbh_src)[i][j][k])
  {
    if ((*nbh)[1][j][k] || (*nbh)[i][1][k] || (*nbh)[i][j][1])
    {
      (*nbh)[i][j][k] = 1;
    }
  }
  i = 0;
  j = 2;
  k = 2;
  if ((*nbh_src)[i][j][k])
  {
    if ((*nbh)[1][j][k] || (*nbh)[i][1][k] || (*nbh)[i][j][1])
    {
      (*nbh)[i][j][k] = 1;
    }
  }
  i = 2;
  j = 0;
  k = 0;
  if ((*nbh_src)[i][j][k])
  {
    if ((*nbh)[1][j][k] || (*nbh)[i][1][k] || (*nbh)[i][j][1])
    {
      (*nbh)[i][j][k] = 1;
    }
  }
  i = 2;
  j = 0;
  k = 2;
  if ((*nbh_src)[i][j][k])
  {
    if ((*nbh)[1][j][k] || (*nbh)[i][1][k] || (*nbh)[i][j][1])
    {
      (*nbh)[i][j][k] = 1;
    }
  }
  i = 2;
  j = 2;
  k = 0;
  if ((*nbh_src)[i][j][k])
  {
    if ((*nbh)[1][j][k] || (*nbh)[i][1][k] || (*nbh)[i][j][1])
    {
      (*nbh)[i][j][k] = 1;
    }
  }
  i = 2;
  j = 2;
  k = 2;
  if ((*nbh_src)[i][j][k])
  {
    if ((*nbh)[1][j][k] || (*nbh)[i][1][k] || (*nbh)[i][j][1])
    {
      (*nbh)[i][j][k] = 1;
    }
  }

  if (nbh_dst)
  {
    for (k = 0; k < 3; k++)
    {
      for (j = 0; j < 3; j++)
      {
        for (i = 0; i < 3; i++)
        {
          (*nbh_dst)[i][j][k] = (*nbh)[i][j][k];
        }
      }
    }
    free(nbh);
    nbh = nbh_dst;
  }

  return nbh;
}

NBH *
N_18_1(NBH * nbh_src, NBH * nbh_dst)
{
  int   i, j, k;
  NBH * nbh;

  if (!nbh_dst)
  {
    nbh = (NBH *)calloc(1, sizeof(NBH));
  }
  else
  {
    nbh = nbh_dst;
  }

  (*nbh)[0][0][0] = TMP;
  (*nbh)[2][0][0] = TMP;
  (*nbh)[0][0][2] = TMP;
  (*nbh)[2][0][2] = TMP;
  (*nbh)[0][2][0] = TMP;
  (*nbh)[2][2][0] = TMP;
  (*nbh)[0][2][2] = TMP;
  (*nbh)[2][2][2] = TMP;
  (*nbh)[1][1][1] = TMP;
  for (k = 0; k < 3; k++)
  {
    for (j = 0; j < 3; j++)
    {
      for (i = 0; i < 3; i++)
      {
        if ((*nbh)[i][j][k] != TMP)
        {
          if ((*nbh_src)[i][j][k])
          {
            (*nbh)[i][j][k] = 1;
          }
          else
          {
            (*nbh)[i][j][k] = 0;
          }
        }
        else
        {
          (*nbh)[i][j][k] = 0;
        }
      }
    }
  }

  return nbh;
}

NBH *
N_18_2(NBH * nbh_src, NBH * nbh_dst)
{
  int   i, j, k;
  NBH * nbh;

  nbh = N_18_1(nbh_src, nullptr);

  i = 0;
  j = 0;
  k = 0;
  if ((*nbh_src)[i][j][k])
  {
    if ((*nbh)[1][j][k] || (*nbh)[i][1][k] || (*nbh)[i][j][1] || (*nbh)[1][1][k] || (*nbh)[1][j][1] || (*nbh)[i][1][1])
    {
      (*nbh)[i][j][k] = 1;
    }
  }
  i = 0;
  j = 0;
  k = 2;
  if ((*nbh_src)[i][j][k])
  {
    if ((*nbh)[1][j][k] || (*nbh)[i][1][k] || (*nbh)[i][j][1] || (*nbh)[1][1][k] || (*nbh)[1][j][1] || (*nbh)[i][1][1])
    {
      (*nbh)[i][j][k] = 1;
    }
  }
  i = 0;
  j = 2;
  k = 0;
  if ((*nbh_src)[i][j][k])
  {
    if ((*nbh)[1][j][k] || (*nbh)[i][1][k] || (*nbh)[i][j][1] || (*nbh)[1][1][k] || (*nbh)[1][j][1] || (*nbh)[i][1][1])
    {
      (*nbh)[i][j][k] = 1;
    }
  }
  i = 0;
  j = 2;
  k = 2;
  if ((*nbh_src)[i][j][k])
  {
    if ((*nbh)[1][j][k] || (*nbh)[i][1][k] || (*nbh)[i][j][1] || (*nbh)[1][1][k] || (*nbh)[1][j][1] || (*nbh)[i][1][1])
    {
      (*nbh)[i][j][k] = 1;
    }
  }
  i = 2;
  j = 0;
  k = 0;
  if ((*nbh_src)[i][j][k])
  {
    if ((*nbh)[1][j][k] || (*nbh)[i][1][k] || (*nbh)[i][j][1] || (*nbh)[1][1][k] || (*nbh)[1][j][1] || (*nbh)[i][1][1])
    {
      (*nbh)[i][j][k] = 1;
    }
  }
  i = 2;
  j = 0;
  k = 2;
  if ((*nbh_src)[i][j][k])
  {
    if ((*nbh)[1][j][k] || (*nbh)[i][1][k] || (*nbh)[i][j][1] || (*nbh)[1][1][k] || (*nbh)[1][j][1] || (*nbh)[i][1][1])
    {
      (*nbh)[i][j][k] = 1;
    }
  }
  i = 2;
  j = 2;
  k = 0;
  if ((*nbh_src)[i][j][k])
  {
    if ((*nbh)[1][j][k] || (*nbh)[i][1][k] || (*nbh)[i][j][1] || (*nbh)[1][1][k] || (*nbh)[1][j][1] || (*nbh)[i][1][1])
    {
      (*nbh)[i][j][k] = 1;
    }
  }
  i = 2;
  j = 2;
  k = 2;
  if ((*nbh_src)[i][j][k])
  {
    if ((*nbh)[1][j][k] || (*nbh)[i][1][k] || (*nbh)[i][j][1] || (*nbh)[1][1][k] || (*nbh)[1][j][1] || (*nbh)[i][1][1])
    {
      (*nbh)[i][j][k] = 1;
    }
  }

  if (nbh_dst)
  {
    for (k = 0; k < 3; k++)
    {
      for (j = 0; j < 3; j++)
      {
        for (i = 0; i < 3; i++)
        {
          (*nbh_dst)[i][j][k] = (*nbh)[i][j][k];
        }
      }
    }
    free(nbh);
    nbh = nbh_dst;
  }

  return nbh;
}

NBH *
N_26_1(NBH * nbh_src, NBH * nbh_dst)
{
  int   i, j, k;
  NBH * nbh;

  if (!nbh_dst)
  {
    nbh = (NBH *)calloc(1, sizeof(NBH));
  }
  else
  {
    nbh = nbh_dst;
  }
  for (k = 0; k < 3; k++)
  {
    for (j = 0; j < 3; j++)
    {
      for (i = 0; i < 3; i++)
      {
        if ((*nbh_src)[i][j][k])
        {
          (*nbh)[i][j][k] = 1;
        }
        else
        {
          (*nbh)[i][j][k] = 0;
        }
      }
    }
  }

  (*nbh)[1][1][1] = 0;

  return nbh;
}

NBH *
Nnk(NBH * nbh_src, NBH * nbh_dst, int connectivity)
{
  NBH * nbh;

  if (!nbh_dst)
  {
    nbh = (NBH *)calloc(1, sizeof(NBH));
  }
  else
  {
    nbh = nbh_dst;
  }

  switch (connectivity)
  {
    case 1: // T6+(x,X)=#C6(N_6_3(x,X))
      N_6_3(nbh_src, nbh);
      break;
    case 2: // T6+(x,X)=#C18(N_18_2)
      N_18_2(nbh_src, nbh);
      break;
    case 3: // T6(x,X)=#C6(N_6_2)
      N_6_2(nbh_src, nbh);
      break;
    case 4: // T26(x,X)=#C6(N_26_1)
      N_26_1(nbh_src, nbh);
      break;
    default: // T6+(x,X)=#C6(N_6_3(x,X))
      N_6_3(nbh_src, nbh);
      break;
  }

  return nbh;
}

/*This function computes the topological number associated with NBH and a certain connectivity.
This function is a "hack": a more elegant/efficient implementation
should be done...
NBH is a binary object (0-1)
connectivity should be equal to 1,2,3, or 4 according to the
convention above*/
int
checkTn(NBH * nbh_src, NBH * nbh_dst, int connectivity)
{
  int i, j, k, a, b, c, ik, jk, kk, ct;
  int con, nvox, label, sum;
  int comp_table[MAX_COMP];
  int min_val;
  int x, y, z;

  NBH * nbh;

  if (!nbh_dst)
  {
    nbh = (NBH *)calloc(1, sizeof(NBH));
  }
  else
  {
    nbh = nbh_dst;
  }

  con = connectivityNumber(connectivity);
  Nnk(nbh_src, nbh, connectivity);

  memset(comp_table, 0, MAX_COMP * sizeof(int));
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      for (k = 0; k < 3; k++)
      {
        if ((*nbh)[i][j][k])
        {
          for (nvox = 0, ik = -1; ik <= 1; ik++)
          {
            a = i + ik;
            if (a < 0 || a >= 3)
            {
              continue;
            }
            for (jk = -1; jk <= 1; jk++)
            {
              b = j + jk;
              if (b < 0 || b >= 3)
              {
                continue;
              }
              for (kk = -1; kk <= 1; kk++)
              {
                sum = abs(ik) + abs(jk) + abs(kk);
                if (sum > con || (!sum))
                {
                  continue;
                }
                c = k + kk;
                if (c < 0 || c >= 3)
                {
                  continue;
                }
                label = (*nbh)[a][b][c];
                if (label > 1)
                {
                  comp_table[label - 1] = 2;
                  nvox++;
                }
              }
            }
          }
          if (!nvox) // find new basin!
          {
            for (ct = 1; comp_table[ct] && ct < MAX_COMP; ct++)
            {
              ;
            }
            (*nbh)[i][j][k] = ct + 1; // label the new basin
            comp_table[ct] = 1;       // note that this number is taken
          }
          else
          {
            min_val = MAX_COMP + 1;
            // merging into the smallest value
            for (ct = 1; ct < MAX_COMP; ct++)
            {
              if (comp_table[ct] == 2)
              {
                min_val = ct;
                break;
              }
            }

            (*nbh)[i][j][k] = min_val + 1;
            comp_table[min_val] = 1;
            // merging of the other neighboring values into the smallest one
            for (ct = min_val + 1; ct < MAX_COMP; ct++)
            {
              if (comp_table[ct] == 2)
              {
                for (x = 0; x < 3; x++)
                {
                  for (y = 0; y < 3; y++)
                  {
                    for (z = 0; z < 3; z++)
                    {
                      if ((*nbh)[x][y][z] == ct + 1)
                      {
                        (*nbh)[x][y][z] = min_val + 1;
                      }
                    }
                  }
                }
                // specify that this basin nbr
                // doesn't exist anymore
                comp_table[ct] = 0;
              }
            }
          }
        }
      }
    }
  }
  for (nvox = 0, ct = 1; ct < MAX_COMP; ct++)
  {
    if (comp_table[ct])
    {
      nvox++;
    }
  }

  if (!nbh_dst)
  {
    free(nbh);
  }

  return nvox;
}

/*This function checks if a point is simple or not
NBH is a binary object (0-1)
connectivity should be equal to 1,2,3, or 4 according to the
convention above*/
int
checkSimple(NBH * nbh, int connectivity)
{
  NBH bnbh, dnbh;

  reverseNBH(nbh, &bnbh);
  if (checkTn(nbh, &dnbh, connectivity) == 1 && checkTn(&bnbh, &dnbh, associatedConnectivity(connectivity)) == 1)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

#endif
