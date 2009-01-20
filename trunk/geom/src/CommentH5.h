// $Id: 
/*
    Copyright (C) 2004 John B. Shumway, Jr.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#ifndef __CommentH5_h_
#define __CommentH5_h_

#include <string>
#include "StructH5.h"

/** The Comment.
    Copies the struct.xml files to a string attribute called comment in struct.h5.
    @author Matthew Harowitz */

class CommentH5 : public StructH5 {
public:
  /** Default Constructor */
  CommentH5();
  /** Constructor */
  CommentH5(const std::string& filename);
  /** Write the comment to a struct.h5 file. */
  void h5Write(const std::string& filename, const int mode=APPEND) const;
  /** The Comment */
  std::string comment;
};
#endif
