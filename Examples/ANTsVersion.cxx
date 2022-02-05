/*=========================================================================
 *  Modified from ANTs reference
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include "ANTsVersion.h"
#include "ANTsVersionConfig.h"

#include <iostream> // std::cout, std::ios
#include <sstream>  // std::ostringstream

namespace
{

std::string
MakeExtendedVersionString()
{
  std::ostringstream v;
  v << "ANTs Version: " << ANTs::Version::VersionString() << std::endl
    << "Compiled: " << ANTs::Version::BuildDate() << std::endl;
  return v.str();
}

static const std::string extendedVersionString = MakeExtendedVersionString();

} // namespace

namespace ANTs
{
unsigned int
Version::MajorVersion()
{
  return ANTS_VERSION_MAJOR;
}
unsigned int
Version::MinorVersion()
{
  return ANTS_VERSION_MINOR;
}
unsigned int
Version::PatchVersion()
{
  return ANTS_VERSION_PATCH;
}
unsigned int
Version::TweakVersion()
{
  return 0;
}
const std::string &
Version::VersionString()
{
  static const std::string v(ANTS_VERSION);
  return v;
}
const std::string &
Version::BuildDate()
{
  static const std::string v(__DATE__ " " __TIME__);
  return v;
}
const std::string &
Version::ExtendedVersionString()
{
  return extendedVersionString;
}
} // namespace ANTs
