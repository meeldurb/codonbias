from Tkinter import *
import mysql.connector
from tkintertable import TableCanvas, TableModel

class oligo_table(Frame):
    def __init__(self):
        Frame.__init__(self)
        self.master.title('oligo')

        self.define_widgets()
        self.grid()

    def define_widgets(self):
        # build table
        root = Tk()
        #attributes = ['oligo_ID', 'oligo_name', 'oligo_type', 'sequence', description, entry_date, creator, update_date, modifier, label5prime, label3prime, labelM1, labelM1position, pathogen_name, target, notes]
        cnx = mysql.connector.connect(user='root', password='root', host='127.0.0.1', database='groupwork')
        cursor = cnx.cursor()
        height = self.count_rows_table()
        width = 14
        query = ("SELECT oligo_ID, oligo_name, oligo_type, sequence, description, entry_date, creator, update_date, modifier, label5prime, label3prime, labelM1, labelM1position, notes FROM Oligo")
        cursor.execute(query)
        for (oligo_ID, oligo_name, oligo_type, sequence, description, entry_date, creator, update_date, modifier, label5prime, label3prime, labelM1, labelM1position, notes) in cursor:
            field_values = [oligo_ID, oligo_name, oligo_type, sequence, description, entry_date, creator, update_date, modifier, label5prime, label3prime, labelM1, labelM1position, notes]
            for i in range(height): #Rows
                for j in range(width): #Columns
                    b = Label(root, text= str(field_values[j]))
                    b.grid(row=i, column=j)
        cursor.close()
        cnx.close()


        



    def make_table(self):
        root = Tk()
        #attributes = ['oligo_ID', 'oligo_name', 'oligo_type', 'sequence', description, entry_date, creator, update_date, modifier, label5prime, label3prime, labelM1, labelM1position, pathogen_name, target, notes]

        cnx = mysql.connector.connect(user='root', password='root', host='127.0.0.1', database='groupwork')
        cursor = cnx.cursor()
        height = self.count_rows_table()
        width = 16
        query = ("SELECT oligo_ID, oligo_name, oligo_type, sequence, description, entry_date, creator, update_date, modifier, label5prime, label3prime, labelM1, labelM1position, notes FROM Oligo")
        cursor.execute(query)
        for (oligo_ID, oligo_name, oligo_type, sequence, description, entry_date, creator, update_date, modifier, label5prime, label3prime, labelM1, labelM1position, notes) in cursor:
            row = [oligo_ID, oligo_name, oligo_type, sequence, description, entry_date, creator, update_date, modifier, label5prime, label3prime, labelM1, labelM1position, notes]
            for i in range(height): #Rows
                for j in range(width): #Columns
                    b = Label(root, text= row[j])
                    b.grid(row=i, column=j)
        cursor.close()
        cnx.close()

    def count_rows_table(self): #count rows in the table in order to know how big the table has to be
        cnx = mysql.connector.connect(user='root', password='root', host='127.0.0.1', database='groupwork')
        cursor = cnx.cursor()
        query = ("SELECT COUNT(oligo_ID) FROM Oligo")
        cursor.execute(query)
        count = self.count_cursor(cursor)
        cursor.close()
        cnx.close()
        return count[0]

    def count_cursor(self, cursor): #retrieves the number stored in cursor, which is needed in count_rows_table()
        for count in cursor:
            return count

if __name__ == "__main__":
    my_gui = oligo_table()
    my_gui.mainloop()
