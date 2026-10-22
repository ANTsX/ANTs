#include "itkantsRegistrationHelper.h"

#include <iostream>
#include <limits>
#include <sstream>

int
main()
{
  using Helper = ants::RegistrationHelper<double, 3>;
  using Image = Helper::ImageType;
  const auto image = Image::New();
  image->SetRegions(Image::SizeType{ { 852, 1182, 1965 } });
  // Geometry alone is sufficient for preflight; do not allocate gigabytes of pixels.
  const auto makeHelper = [&](Helper::XfrmMethod method, bool addImageMetric = true) {
    auto helper = Helper::New();
    std::vector<unsigned int> mesh(3, 4);
    switch (method)
    {
      case Helper::Rigid:
        helper->AddRigidTransform(0.1);
        break;
      case Helper::Affine:
        helper->AddAffineTransform(0.1);
        break;
      case Helper::BSpline:
        helper->AddBSplineTransform(0.1, mesh);
        break;
      case Helper::GaussianDisplacementField:
        helper->AddGaussianDisplacementFieldTransform(0.1, 3, 0);
        break;
      case Helper::BSplineDisplacementField:
        helper->AddBSplineDisplacementFieldTransform(0.1, mesh, mesh, 3);
        break;
      case Helper::SyN:
        helper->AddSyNTransform(0.1, 3, 0);
        break;
      case Helper::BSplineSyN:
        helper->AddBSplineSyNTransform(0.1, mesh, mesh, 3);
        break;
      case Helper::TimeVaryingVelocityField:
        helper->AddTimeVaryingVelocityFieldTransform(0.1, 4, 3, 0, 0, 0);
        break;
      case Helper::TimeVaryingBSplineVelocityField:
        helper->AddTimeVaryingBSplineVelocityFieldTransform(0.1, mesh, 4, 3);
        break;
      case Helper::Exponential:
        helper->AddExponentialTransform(0.1, 3, 0, 4);
        break;
      case Helper::BSplineExponential:
        helper->AddBSplineExponentialTransform(0.1, mesh, mesh, 4, 3);
        break;
      default:
        break;
    }
    if (addImageMetric)
    {
      helper->AddMetric(Helper::CC, image, image, 0, 1, Helper::random, 32, 4, false, 0.01);
    }
    helper->SetIterations({ { 10, 10 } });
    helper->SetShrinkFactors({ { 2, 1 } });
    helper->SetSmoothingSigmas({ { 1, 0 } });
    helper->SetSmoothingSigmasAreInPhysicalUnits({ false });
    return helper;
  };
  int failures = 0;
  const auto expect = [&](bool condition, const char * message) {
    if (!condition)
    {
      std::cerr << message << '\n';
      ++failures;
    }
  };
  const bool limited = std::numeric_limits<Helper::ObjectMetricType::NumberOfParametersType>::max() / 3 <
                       852ULL * 1182 * 1965;
  for (const auto method : { Helper::GaussianDisplacementField, Helper::BSplineDisplacementField, Helper::SyN,
                            Helper::BSplineSyN, Helper::TimeVaryingVelocityField,
                            Helper::TimeVaryingBSplineVelocityField, Helper::Exponential, Helper::BSplineExponential })
  {
    auto helper = makeHelper(method);
    std::ostringstream log;
    helper->SetLogStream(log);
    auto * previousErrorBuffer = std::cerr.rdbuf(log.rdbuf());
    const auto result = helper->ValidateParameters();
    std::cerr.rdbuf(previousErrorBuffer);
    expect(result == (limited ? EXIT_FAILURE : EXIT_SUCCESS), "Dense-field limit not enforced");
    if (limited)
    {
      expect(log.str().find("Stage 0") != std::string::npos && log.str().find("level 2") != std::string::npos,
             "Error does not identify the offending stage and level");
    }
    helper->SetIterations({ { 10, 0 } });
    expect(helper->ValidateParameters() == EXIT_SUCCESS, "Zero-iteration full-resolution level was rejected");
    helper->SetIterations({ { 10, 10 } });
    helper->SetShrinkFactors({ { 4, 2 } });
    expect(helper->ValidateParameters() == EXIT_SUCCESS, "Safe coarser levels were rejected");
  }
  for (const auto method : { Helper::Rigid, Helper::Affine, Helper::BSpline })
  {
    expect(makeHelper(method)->ValidateParameters() == EXIT_SUCCESS, "Non-dense transform was rejected");
  }

  auto masked = makeHelper(Helper::SyN);
  auto mask = Helper::MaskImageType::New();
  mask->SetRegions(Helper::MaskImageType::SizeType{ { 1, 1, 1 } });
  masked->AddFixedImageMask(mask);
  expect(masked->ValidateParameters() == (limited ? EXIT_FAILURE : EXIT_SUCCESS),
         "An image metric mask incorrectly reduced the parameter count");

  // A single point-set metric takes a point-sized derivative path in the two
  // B-spline methods, even when the mask defining the virtual domain is large.
  auto largeMask = Helper::MaskImageType::New();
  largeMask->SetRegions(image->GetLargestPossibleRegion());
  const auto points = Helper::LabeledPointSetType::New();
  Helper::LabeledPointSetType::PointType point;
  point.Fill(1.0);
  points->SetPoint(0, point);
  points->SetPointData(0, 1);
  const auto addPointSetMetric = [&](Helper * helper) {
    helper->AddMetric(Helper::ICP, nullptr, nullptr, points, points, nullptr, nullptr,
                      0, 1, Helper::none, 32, 4, false, false, 1.0, 50, 1.1, false, 1.0,
                      std::sqrt(5.0), std::sqrt(5.0));
  };
  for (const auto method : { Helper::BSplineSyN, Helper::TimeVaryingBSplineVelocityField })
  {
    auto pointOnly = makeHelper(method, false);
    pointOnly->AddFixedImageMask(largeMask);
    addPointSetMetric(pointOnly);
    expect(pointOnly->ValidateParameters() == EXIT_SUCCESS,
           "Single point-set B-spline stage was incorrectly rejected by the voxel limit");

    addPointSetMetric(pointOnly);
    expect(pointOnly->ValidateParameters() == (limited ? EXIT_FAILURE : EXIT_SUCCESS),
           "Multiple point-set metrics incorrectly bypassed the voxel limit");

    auto mixed = makeHelper(method);
    mixed->AddFixedImageMask(largeMask);
    addPointSetMetric(mixed);
    expect(mixed->ValidateParameters() == (limited ? EXIT_FAILURE : EXIT_SUCCESS),
           "Mixed image and point-set metrics incorrectly bypassed the voxel limit");
  }
  for (const auto method : { Helper::SyN, Helper::TimeVaryingVelocityField })
  {
    auto pointOnly = makeHelper(method, false);
    pointOnly->AddFixedImageMask(largeMask);
    addPointSetMetric(pointOnly);
    expect(pointOnly->ValidateParameters() == (limited ? EXIT_FAILURE : EXIT_SUCCESS),
           "Gaussian point-set registration incorrectly bypassed the voxel limit");
  }

  if (limited)
  {
    const auto maximumVoxels = std::numeric_limits<Helper::ObjectMetricType::NumberOfParametersType>::max() / 3;
    image->SetRegions(Image::SizeType{ { maximumVoxels, 1, 1 } });
    expect(makeHelper(Helper::SyN)->ValidateParameters() == EXIT_SUCCESS, "Maximum valid count was rejected");
    image->SetRegions(Image::SizeType{ { maximumVoxels + 1ULL, 1, 1 } });
    expect(makeHelper(Helper::SyN)->ValidateParameters() == EXIT_FAILURE, "First invalid count was accepted");
  }

  // Anisotropic spacing means scalar factor 2 need not shrink all three axes.
  image->SetSpacing(Image::SpacingType{ { 1, 10, 10 } });
  image->SetRegions(Image::SizeType{ { 2000, 2000, 1000 } });
  auto anisotropic = makeHelper(Helper::SyN);
  anisotropic->SetIterations({ { 10, 0 } });
  expect(anisotropic->ValidateParameters() == (limited ? EXIT_FAILURE : EXIT_SUCCESS),
         "Anisotropic shrink factors were not respected");

  // Check a later stage before any earlier optimization can start, with reversed metric IDs.
  image->SetSpacing(Image::SpacingType{ { 1, 1, 1 } });
  image->SetRegions(Image::SizeType{ { 852, 1182, 1965 } });
  auto later = makeHelper(Helper::Affine);
  later->AddSyNTransform(0.1, 3, 0);
  later->AddMetric(Helper::CC, image, image, 1, 1, Helper::none, 32, 4, false, 1);
  later->SetIterations({ { 10 }, { 10 } });
  later->SetShrinkFactors({ { 1 }, { 1 } });
  later->SetSmoothingSigmas({ { 0 }, { 0 } });
  later->SetSmoothingSigmasAreInPhysicalUnits({ false, false });
  expect(later->ValidateParameters() == (limited ? EXIT_FAILURE : EXIT_SUCCESS), "Later dense-field stage not checked");

  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
