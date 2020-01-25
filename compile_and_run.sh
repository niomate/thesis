echo "Compiling..."
cmake --build build
echo "Finished compiling, running test..."
./build/test &> output.log
echo "Finished testing, view log ? [Y/n]"
read decision
if [ $decision = "y" ]; then
   less output.log
fi
#echo "Run chain visualization? [Y/n]"
#read decision
#if [ $decision = "y" ]; then
   #python scripts/visualize_chains.py
#fi
