#ifndef ANTSCOUT_HXX
#define ANTSCOUT_HXX

#include <cassert>
// #include "boost/iostreams/stream.hpp"
// #include "boost/iostreams/concepts.hpp"
#include <iostream>

namespace ants
{

/*
class ants_Sink : public boost::iostreams::sink
{
public:
  explicit ants_Sink() : os_( NULL )
  {
  }

  void set_stream( std::ostream* os )
  {
    // assert failed means output streams is being being changed within program
    assert( os_ == NULL );
    os_ = os;
  }

  std::streamsize write( const char* buffer, std::streamsize num_chars )
  {
    if( os_ != NULL )
      {
      os_->write( buffer, num_chars );
      os_->flush();
      }
    return num_chars;
  }

private:
  // user provided output stream; defaults to NULL
  std::ostream* os_;
};
*/
} // namespace ants

#endif // ANTSCOUT_HXX
