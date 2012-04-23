/**
 * There were several functions that had been copied and
 * pasted over and over again in this library.  This
 * header files is contains a common definitoin for
 * those file.
 * \author Hans J. Johnson
 */
#ifndef __antsUtilities_h__
#define __antsUtilities_h__

#include "antscout.hxx"
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include <errno.h>

// We need to ensure that only one of these exists!
namespace ants
{
extern boost::iostreams::stream<ants_Sink> antscout;
}

// ##########################################################################
// TODO: KENT:  This block feels like it could be better encapsulated as a c++ class
//
typedef enum
  {
  INVALID_FILE = 1,
  AFFINE_FILE,
  DEFORMATION_FILE,
  IMAGE_AFFINE_HEADER,
  IDENTITY_TRANSFORM
  } TRAN_FILE_TYPE;

// TODO: This should be a class.
typedef struct
  {
  //    char *filename;
  std::string filename;
  TRAN_FILE_TYPE file_type;
  bool do_affine_inv;
  //    void SetValue(char *filename, TRAN_FILE_TYPE file_type, bool do_affine_inv){
  //        this.filename = filename;
  //        this.file_type = file_type;
  //        this.do_affine_inv = do_affine_inv;
  //    };
  double weight;   // for average
  } TRAN_OPT;

typedef std::vector<TRAN_OPT> TRAN_OPT_QUEUE;

typedef struct
  {
  bool physical_units;
  std::vector<double> sigma;
  } MLINTERP_OPT;

typedef struct
  {
  bool use_NN_interpolator;
  bool use_MultiLabel_interpolator;
  bool use_BSpline_interpolator;
  bool use_TightestBoundingBox;
  char * reference_image_filename;
  bool use_RotationHeader;

  MLINTERP_OPT opt_ML;
  } MISC_OPT;

extern TRAN_FILE_TYPE CheckFileType(const char * const str);

extern TRAN_FILE_TYPE CheckFileType(const std::string str);

extern void SetAffineInvFlag(TRAN_OPT & opt, bool & set_current_affine_inv);

extern void DisplayOptQueue(const TRAN_OPT_QUEUE & opt_queue);

extern void DisplayOpt(const TRAN_OPT & opt);

// ##########################################################################

extern bool get_a_double_number(const char * const str, double & v);

// TODO: KENT:  These two functions have cross-platform-equivalent versions from kwSys and could be replaced.
extern void FilePartsWithgz(const std::string & filename, std::string & path, std::string & name, std::string & ext);

extern bool CheckFileExistence(const char * const str);

#endif // __antsUtilities_h__
