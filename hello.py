import os
import boto3
 
print("File upload starting")
 
s3 = boto3.client("s3")
bucket_name = "binf5506-hoang"
 
s3.upload_file("hello.py", bucket_name, "hello.py")
 
print("File upload completed")