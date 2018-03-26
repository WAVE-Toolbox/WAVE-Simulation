#!/bin/bash

# Change to src/ folder
cd ../

# Format *.cpp files
<<<<<<< HEAD
find $(pwd) -name '*.cpp' -exec clang-format-mp-3.9 -style=file -i {} \;
=======
find $(pwd) -name '*.cpp' -exec clang-format -style=file -i {} \;
>>>>>>> origin/develop

# Format *.hpp files
find $(pwd) -name '*.hpp' -exec clang-format-mp-3.9 -style=file -i {} \;
