import psycopg2
import pprint as pp
# load the psycopg extras module
import psycopg2.extras
from psycopg2.extras import RealDictCursor

# Connect to an existing database
conn = psycopg2.connect(host="ec2-52-14-42-130.us-east-2.compute.amazonaws.com", port="5432", dbname="postgres", user="postgres", password="postgres")

# Open a cursor to perform database operations
cur = conn.cursor(cursor_factory=RealDictCursor)

# Execute a command: this creates a new table
cur.execute("""SELECT * FROM prodtest""""")


# Pass data to fill a query placeholders and let Psycopg perform
# the correct conversion (no more SQL injections!)
records = cur.fetchall()

#print (records)

for row in records:
    for col in row:
        print (col, row[col])
     
    
# Query the database and obtain data as Python objects

column_names = []
data_rows = []

column_names = [desc[0] for desc in cur.description]

for row in cur:
    data_rows.append(row)

print("Column names: {}\n".format(column_names))

#print cur
# Make the changes to the database persistent
#conn.commit()

# Close communication with the database
cur.close()
conn.close()


