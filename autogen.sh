#!/bin/bash
libtoolize && aclocal && autoconf && automake --add-missing
echo "Now run ./configure to configure the program"
echo "You can view relevant options by running ./configure -h"
