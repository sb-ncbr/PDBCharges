import sys

def replace_string(file_path, search_string, replacement_string):
    try:
        # Read the file
        with open(file_path, 'r') as file:
            content = file.read()

        # Count occurrences of the search string
        occurrences = content.count(search_string)

        if occurrences == 0:
            print(f"Error: String '{search_string}' not found in file '{file_path}'.")
            sys.exit(1)
        elif occurrences > 1:
            print(f"Error: String '{search_string}' occurs multiple times ({occurrences}) in file '{file_path}'.")
            sys.exit(1)

        # Replace the string
        new_content = content.replace(search_string, replacement_string, 1)

        # Write the modified content back to the file
        with open(file_path, 'w') as file:
            file.write(new_content)

        print(f"String '{search_string}' successfully replaced with '{replacement_string}'.")

    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        sys.exit(1)

# Ensure all arguments are provided
if len(sys.argv) != 2:
    print("Usage: python edit_moleculekit.py <path_to_preparation.py>")
    sys.exit(1)

file_path = sys.argv[1]

# Replace strings in the file preparation.py of Moleculekit version 1.9.15
replace_string(file_path,
    '("charge", "formalcharge"),',
    '("charge", "formalcharge"),\n        ("ffcharge", "charge"),'
)
replace_string(file_path,' opt=True,', ' opt=False,')

