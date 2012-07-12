# \author Hans J. Johnson
# This file will help store new data
# in a publicly available web location
# for the purpose of allowing larger
# testing data to be separate from
# the source tree.
#
# NOTE:  This is not a global solution
#        for any developer.  Only members
#        of the Iowa PINC lab have access
#        to generate data in this way.
#
#
import hashlib
import os
import argparse
import sys
import shutil

def md5_for_file(f, block_size=2**20):
  """Generate a hash key from a file"""
  md5 = hashlib.md5()
  while True:
      data = f.read(block_size)
      if not data:
          break
      md5.update(data)
  return md5.hexdigest()


if __name__ == '__main__':
  defaultCMakeRepository='/iplweb/html/users/brainstestdata/ctestdata'

  parser = argparse.ArgumentParser(description="Program to move local test data to external repository")
  parser.add_argument('--src',dest='sourceFile',help='The source file to be moved.')
  parser.add_argument('--dest',dest='destPath',default=defaultCMakeRepository,help='The external repository location')

  args = parser.parse_args()

  fileName = args.sourceFile
  # Do not upload file hashes of file hashs.  To prevent infinite recursion when using sloppy glob expressions!
  if fileName[-4:] == ".md5":
      print("Skipping md5 file:  {0}".format(fileName))
      exit(0)
  print("Preparing: {0}".format(fileName))
  localPublicPath = args.destPath
  algo = 'MD5'
  f = open(fileName)
  value = md5_for_file(f)
  f.close()
  #print('MD5 value={0}'.format(value))
  #print('MD5 value={0}'.format('59871b1b19b16a3c04d752f54bbf8bfd'))
  source = fileName
  destPath=localPublicPath+'/'+algo
  if not os.path.exists(destPath):
      os.mkdir(destPath)
  dest = destPath +'/'+value
  print("mv -f {0} {1}".format(source,dest))
  md5FileName=fileName+'.md5'
  f = open(md5FileName,'w')
  f.write(value)
  f.close()
  finalDestination=destPath+'/'+os.path.basename(md5FileName)
  if os.path.exists(finalDestination):
    print("Destination file already exists: SKIPPING {0}".format(finalDestination))
  else:
    shutil.copyfile(md5FileName,destPath+'/'+os.path.basename(md5FileName))
    shutil.copy(source,dest)
    os.unlink(source)
  ## if prepareing data remotely, echo a helpful rsync command needed to push from remote destination to IPL via rsync
  if args.destPath != defaultCMakeRepository:
    print('rsync -av {0}/ neuron.psychiatry.uiowa.edu:{1}/'.format(args.destPath,defaultCMakeRepository))
