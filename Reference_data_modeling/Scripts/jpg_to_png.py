import sys
import os
from PIL import Image

def get_file_paths(top_directory, extension):
    imageFiles = []

    for root, _, files in os.walk(top_directory):
        for file in files:
            if file.endswith(extension):
                imageFiles.append(os.path.join(root, file))
                
    return imageFiles

top_directory = r'/media/mc/2TB/Dropbox/UUV Project/Pictures/Surface Area Photos'

# get all the jpg file paths
jpg_files = get_file_paths(top_directory, '.jpg')

for file in jpg_files:
    root_path = os.path.split(file)[0]
    file_name = os.path.split(file)[1][:-4]
    try:
        img = Image.open(file)
        img.save(os.path.join(root_path, file_name + '.png'))
        os.remove(file)
    except FileNotFoundError:
        print('Could not convert %s' % file)

print("File conversion completed!")