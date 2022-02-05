/*=========================================================================
 *
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
#ifndef ANTsVersion_h
#define ANTsVersion_h

// NOTE: Windows magic needed here
#define ANTsCommon_EXPORT

#include <string>

namespace ANTs
{

/** \class Version
 * \brief Version info for ANTs
 */
class ANTsCommon_EXPORT Version
{
public:
  static unsigned int
  MajorVersion();
  static unsigned int
  MinorVersion();
  static unsigned int
  PatchVersion();
  static unsigned int
  TweakVersion();
  static const std::string &
  VersionString();
  static const std::string &
  BuildDate();

  static const std::string &
  ExtendedVersionString();
  std::string
  ToString()
  {
    return Version::ExtendedVersionString();
  }
};
} // namespace ANTs

#endif
