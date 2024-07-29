def log_msg(msg, filename):
    with open(filename, 'a') as file:
        file.write(msg)
    print(msg)