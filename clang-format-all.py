import subprocess
import os

files = [f for f in os.listdir(".") if os.path.isfile(f)]
# print(files)
files_filtered = [f for f in files if ".cpp" in f]
# print(files_filtered)
# clang-format -i Virtual_Tree.cpp
for f in files_filtered:
    print("Formatting", f)
    try:
        subprocess.run(["clang-format", "-i", f])
        print("Successfully Formated", f)
    except Exception as e:
        print("Failed Formated", f)
        print(e)
