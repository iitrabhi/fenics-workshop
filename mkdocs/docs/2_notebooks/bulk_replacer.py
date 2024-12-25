import os

# Base directory containing day folders
base_dir = os.getcwd() 

# Iterate through day folders
for day in [f"day-{i}" for i in range(1, 6)]:
    day_path = os.path.join(base_dir, day)
    
    for nature in ["exercises", "tutorials"]:
        nature_path = os.path.join(day_path, nature)
        
        if os.path.exists(nature_path):
            # Iterate through Markdown files
            for file_name in os.listdir(nature_path):
                if file_name.endswith(".md"):
                    file_path = os.path.join(nature_path, file_name)
                    file_stem = os.path.splitext(file_name)[0]  # Filename without extension
                    
                    # Prepare the line to add
                    new_line = f"The accompanying Jupyter notebook can be obtained here [{file_stem}](../../../../../src/{day}/{nature}/{file_stem}.ipynb)\n\n"
                    
                    # Read the current content of the file
                    with open(file_path, "r", encoding="utf-8") as file:
                        content = file.read()
                    
                    # Write the new line followed by the existing content
                    with open(file_path, "w", encoding="utf-8") as file:
                        file.write(new_line + content)
                    
                    print(f"Updated: {file_path}")
