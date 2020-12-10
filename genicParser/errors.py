
def path_invalid(parent_path, operation):
    return f"INVALID Path for {operation}\n" \
           f"{operation} attempt to navigate to a directory or file at {parent_path} but it does not exist"