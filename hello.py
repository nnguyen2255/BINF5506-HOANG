import os
import boto3
 
s3 = boto3.client("s3")
bucket_name = "sohail-binf55062"
 
s3.upload_file("hello.py", bucket_name, "hello.py")
 
print("File upload completed")