mkdir out

python3 timings.py --threads 1  --example 1 --order 5 --mesh_size 0.1 -def 0
python3 timings.py --threads 2  --example 1 --order 5 --mesh_size 0.1 -def 0
python3 timings.py --threads 4  --example 1 --order 5 --mesh_size 0.1 -def 0
python3 timings.py --threads 12 --example 1 --order 5 --mesh_size 0.1 -def 0

python3 timings.py --threads 1  --example 2 --order 4 --mesh_size 0.125 -sm 1 -def 0
python3 timings.py --threads 2  --example 2 --order 4 --mesh_size 0.125 -sm 1 -def 0
python3 timings.py --threads 4  --example 2 --order 4 --mesh_size 0.125 -sm 1 -def 0
python3 timings.py --threads 12 --example 2 --order 4 --mesh_size 0.125 -sm 1 -def 0