import tkinter

class myTime:
    def __init__(self, str_time):
        if str_time.count(':') == 1:
            str_time = "00:" + str_time
        print(str_time)
        time = str_time.split(':')
        self.abs_seconds = int(time[2]) + int(time[1]) * 60
        + int(time[0]) * 3600
        
    def hours(self):
        return self.abs_seconds // 3600 
        
    def minutes(self):
        return (self.abs_seconds % 3600) / 60
    
    def seconds(self):
        return self.abs_seconds % 60
        
    def __sub__(self, other):
        return self.abs_seconds - other.seconds
    
    def __le__(self, other):
        return self.abs_seconds <= other.seconds
    
    def __str__(self):
        return str(self.seconds) + ":" + str(self.minutes) + ":" + str(self.seconds)
        

def tb_input(tb_text, tb2):
    pass

def set_text():
    get_breaks()
    new_text = get_time_tags()
    text_box1.delete("1.0", tkinter.END)
    text_box1.insert("1.0", new_text)
    
def get_time_tags():
    strings = text_box1.get("1.0", 'end-1c').split("\n")
    result = []
    tmp = []
    for string in strings:
        tmp = string.split()
        result.append([myTime(tmp[0]), ' '.join(tmp[1:])])
    return result

def get_breaks():
    strings = text_box2.get("1.0", 'end-1c').split("\n")
    result = []
    for string in strings:
        tmp = string.split()
        for i in range(len(tmp)):
            result.append([myTime(tmp[0]), myTime(tmp[1])])
    return result
        
def change_time_tags():
    time_tags = get_time_tags()
    breaks = get_breaks
    cur_break = 0
    addition = myTime("00:00:00")
    for time_tag in time_tags:
        if breaks[cur_break][0] <= time_tag[0]:
            addition += breaks[cur_breaks][1] - breaks[cur_breaks][0]
            cur_break += 1
        

root = tkinter.Tk()
text_box1 = tkinter.Text(root, width=34, height=15, bg="#FFFFE0")
text_box2 = tkinter.Text(root, width=34, height=5, bg="#F0F8FF")
label1 = tkinter.Label(root, text='time tags:', font="Helvetica 10")
label2 = tkinter.Label(root, text='breaks:', font="Helvetica 10")
button1 = tkinter.Button(root, width=33, text='run', font="Helvetica 10 italic", 
                         command=lambda:set_text())
label1.pack()
text_box1.pack()
label2.pack()
text_box2.pack()
button1.pack()
root.mainloop()