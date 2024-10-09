CHARS_LIST=("@AB" "@AB" "@ABC" "@AB" "@A" "@A" "@" "@AB" "01" "01" "01" "01" "01" "01" "01" "01" "01" "01" "01" "01")  # Example: 1st index uses 'abc', 2nd uses 'xyz', 3rd uses '123'
LENGTH=${#CHARS_LIST[@]}       
PYTHON_SCRIPT="./processing.py"  # Python script to run

# Function to generate permutations with different possible characters for each index
perms() {
    local chars_list=("$@")
    local length=${#chars_list[@]}

    python_chars_list="["
    for chars in "${chars_list[@]}"; do
        python_chars_list+="\"$chars\", "  
    done

    python_chars_list="${python_chars_list%, }]"
    
    python3 -c "
import itertools

chars_list = [list(chars) for chars in $python_chars_list]  # Convert each string into a list of characters

# Generate cartesian product of the character sets
for p in itertools.product(*chars_list):
    print(''.join(p))
    "
}

run_in_parallel() {
    local chars_list=("$@")
    local python_script=$PYTHON_SCRIPT

    # Generate permutations and run the Python script for each in parallel
    for perm in $(perms "${chars_list[@]}"); do
        python3 "$python_script" "$perm" &
    done

    wait
}

run_in_parallel "${CHARS_LIST[@]}"