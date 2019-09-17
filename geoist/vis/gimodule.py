#Load graphical tool kit
#simple GUI framework for geoist

from tkinter import *
import os

########################################
##############  Constants  #############
########################################

maxrow    = 200
maxwidth  = 80
textwidth = 25

color_seperator  = "#5E84A9"
color_mainline   = "#3F5871"
color_button     = "#8EC7FF" #D8DADF"

pythonpath   = os.path.join(os.path.dirname(__file__))
getcwd       = os.getcwd()

####################################
###############  Functions  ########
####################################

def seperator_line(Iframe, Ilineposition, Itext=None):
    line = Label(Iframe, anchor=W, relief=RAISED,
    width=maxwidth, text=Itext, fg="white", bg=color_seperator)
    line.grid(row=Ilineposition, column=0, columnspan=5, sticky=W)



def runprogram(Icommand, Ijobfile, Ioutput):

    if os.name == "posix":
        os.system(pythonpath+"/bin/"+Icommand+"< "+Ijobfile+" | tee "+Ioutput)
        print("##### Program completed #####")
    elif os.name == "nt":
        import subprocess

        outputfile = open(Ioutput, 'w')
        process = subprocess.Popen([os.path.join(pythonpath, "bin", Icommand + '.exe')],
        stdin=open(Ijobfile, "rb"),
        stdout=subprocess.PIPE)

        print("#### Program running... #####")
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() != None: break
            outputfile.writelines(output)
            print(output),
        outputfile.close
        print("##### Program completed #####")

######################################
############   Classes   #############
######################################

class statusbar:
    def config(self,Itext,Itype):
        if Itype == "normal":
            self.statusbar.config(text=Itext, fg="dark green")
        elif Itype == "warning":
            self.statusbar.config(text=Itext, fg="red")
        else:
            self.statusbar.config(text=Itext, fg=Itype)    
    
    def __init__(self, Iframe, Ilineposition):
        self.frame = Iframe
        self.lineposition = Ilineposition
    
        self.statusbar = Label(Iframe, bd=1, relief=SUNKEN, anchor=W, width=maxwidth, text="")
        self.statusbar.grid(row=Ilineposition, column=0, columnspan=4, sticky=W)



class mainline:
    def __init__(self, Iframe, Iprogramname, Iversion):
        mainLabel = Label(Iframe, anchor=W, relief=RAISED,
        text=Iprogramname,
        justify="left", fg="white", bg=color_mainline, width=maxwidth)
        mainLabel.grid(row=0, column=0, columnspan=4, sticky=W)



class liste:
    def __init__(self, Iframe, Ilineposition, Iquestion, Ichoices):
        size = len(Ichoices)
        self.list = Listbox(Iframe, height=size, selectmode="SINGLE",
        takefocus=0, width = textwidth)
        for item in Ichoices:
            self.list.insert(END, item)
        self.list.selection_set(0)
        self.listLabel = Label(Iframe,text=Iquestion)
        self.listLabel.grid(row=Ilineposition, column=0, columnspan=2, sticky=W)
        self.list.grid(row=Ilineposition, column=1, columnspan=2, sticky=W)



class menubar:
    def showhelp(self):
        helpWindows = SimpleDialog(self.root,
        text=self.descrip, buttons=["OK"], default=0, title="Help").go()

    def __init__(self, Iframe=None, Iquit=None, Iwrite=None,
        Irun=None, Ihelp=None, Iroot=None):
        if Iroot != None:
            self.root = Iroot
            self.descrip = Ihelp

        #Initiates the menubars
        menubar = Menu(Iframe)

        # create a pulldown menu, and add it to the menu bar
        filemenu = Menu(menubar, tearoff=0)
        filemenu.add_command(label="Write settings", command=Iwrite)
        filemenu.add_command(label="Run program", command=Irun)
        filemenu.add_separator()
        filemenu.add_command(label="Quit", command=Iquit)
        menubar.add_cascade(label="File", menu=filemenu)

        helpmenu = Menu(menubar, tearoff=0)
        helpmenu.add_command(label="About", command=self.showhelp)
        menubar.add_cascade(label="Help", menu=helpmenu)

        #Draws the menu in the program
        Iroot.config(menu=menubar, relief=RAISED, borderwidth=2)
        #flat, groove, raised, ridge, solid, or sunken


class quitwriterun:
    def showhelp(self):
        helpWindows = SimpleDialog(self.root,
        text=self.descrip, buttons=["OK"], default=0, title="Help").go()

    def __init__(self, Iframe, Ilineposition, Iquit,
        Iwrite, Irun, Ihelp=None, Iroot=None):
        if Iroot != None:
            self.root = Iroot
            self.descrip = Ihelp

        button = Button(Iframe, text="Quit", fg="black", command=Iquit)
        button.grid(row=Ilineposition, column=0, sticky=W)

        writetofileButton = Button(Iframe, fg="black",
        text="Write settings", command=Iwrite)
        writetofileButton.grid(row=Ilineposition, column=1, sticky=W)

        rungeolcolButton = Button(Iframe, fg="black",
        text="Run program", command=Irun)
        rungeolcolButton.grid(row=Ilineposition, column=2, sticky=W)

        if Ihelp != None:
            help = Button(Iframe, text="Help", command=self.showhelp)
            help.grid(row=Ilineposition, column=3, sticky=E)



class NewRadioButton:
    def enable(self):
        self.on.configure( state='normal')
        self.off.configure(state='normal' )
        #self.off.configure(state='normal', selectcolor="red")

    def disable(self):
        self.on.configure( state='disabled')
        self.off.configure(state='disabled')
        #self.off.configure(state='disabled', selectcolor="grey")

    def draw(self, Irow):
        self.label.grid(row=Irow, column=0, sticky=W)
        self.on.grid(row=Irow, column=1, sticky=W)
        self.off.grid(row=Irow, column=2, sticky=W)
        if self.descrip != "nohelp":
            self.help.grid(row=Irow, column=4, sticky=W)

    def showhelp(self):
        helpWindows = SimpleDialog(self.root, text=self.descrip,
        buttons=["OK"], default=0, title="Help").go()

    def __init__(self, Iframe, Iroot, Ilabel, Iselected, IcommandYES, IcommandNO, Idescrip):
        self.varname = Ilabel
        self.state = StringVar()
        self.descrip = Idescrip
        self.root = Iroot

        if IcommandYES == "nocommand":
            self.on = Radiobutton(Iframe, text="Yes", value="t",
                    width=3, selectcolor=color_button, indicatoron=0, variable=self.state )
            self.off = Radiobutton(Iframe, text="No", value="f",
                    width=3, selectcolor=color_button, indicatoron=0, variable=self.state )
        else:
            self.on = Radiobutton(Iframe, text="Yes", value="t",
                    width=3, selectcolor=color_button, indicatoron=0, variable=self.state,
                    command=IcommandYES)
            self.off = Radiobutton(Iframe, text="No", value="f",
                    width=3, selectcolor=color_button, indicatoron=0, variable=self.state,
                    command=IcommandNO)

        self.label = Label(Iframe, text=Ilabel)

        if Idescrip != "nohelp":
            self.help = Button(Iframe, text="?", command=self.showhelp)

        if Iselected == "select":
            self.on.select()
            self.selected  = "t"
        else:
            self.off.select()
            self.selected = "f"

    def invokeyes(self):
        self.on.invoke()

    def invokeno(self):
        self.off.invoke()


class LauncherButton:
    # This is used in the launcher.py to start the programs.

    def run(self, Icommand):
        if os.name == "posix":
            os.system(pythonpath+"/"+Icommand+"&")
        elif os.name == "nt":
            os.system(Icommand)
#	    import subprocess
#	    print os.path.join(pythonpath, Icommand)
#            p = subprocess.Popen([os.path.join(pythonpath, Icommand)],
#                                  stdout=subprocess.PIPE,shell=True)

    def __init__(self, Iframe, Irow, Iwidth, Icolumn, Iquestion,
        Icommand, Iexplanation):
        # Icolumn=0 : First row
        # Icolumn=1 : Second row
        self.button = Button(Iframe, text=Iquestion,
        width=Iwidth, command=Icommand)
        self.label = Label(Iframe, text=Iexplanation)
        self.button.grid( row=Irow, column=0+2*int(Icolumn), sticky=W)
        self.label.grid( row=Irow, column=1+2*int(Icolumn), sticky=W)



class NewEntry:
    def enable(self):
        self.textentry.configure( state='normal')

    def disable(self):
        self.textentry.configure( state='disabled')

    def showhelp(self):
        helpWindows = SimpleDialog(self.root,
        text=self.descrip, buttons=["OK"], default=0, title="Help").go()

    def __init__(self, Iframe, Iroot, Irow,
        Iwidth, Istate, Ilabel, IstartText, Idescrip):
        self.text = StringVar()
        self.descrip = Idescrip
        self.root = Iroot
        self.textentry = Entry(Iframe, width=Iwidth, state=Istate, textvariable=self.text)
        self.textentry.insert(0, IstartText)
        self.label = Label(Iframe, text=Ilabel)
        self.label.grid(row=Irow, column=0, sticky=W)
        self.textentry.grid(row=Irow, column=1, sticky=W, columnspan=2)

        if Idescrip != "nohelp":
            self.help = Button(Iframe, text="?", command=self.showhelp)
            self.help.grid(row=Irow, column=4, sticky=E)



class FileSelector:
    def enable(self):
        self.textentry.configure( state='normal')

    def disable(self):
        self.textentry.configure( state='disabled')

    def showhelp(self):
        helpWindows = SimpleDialog(self.root,
        text=self.descrip, buttons=["OK"], default=0, title="Help").go()

    def showfilebrowser(self):
        import tkFileDialog
        file = tkFileDialog.askopenfile()
        try:
            #Only do something if a file has been selected.
            if file.name:
                self.textentry.delete(0,255)
                self.textentry.insert(0, file.name)
        except:
            #If no file has been selected do nothing.
            pass

    def __init__(self, Iframe, Iroot, Irow,
        Iwidth, Istate, Ilabel, IstartText, Idescrip):
        self.text = StringVar()
        self.descrip = Idescrip
        self.root = Iroot
        self.textentry = Entry(Iframe, width=Iwidth, state=Istate, textvariable=self.text)
        self.textentry.insert(0, IstartText)
        self.label = Label(Iframe, text=Ilabel)
        self.label.grid(row=Irow, column=0, sticky=W)
        self.textentry.grid(row=Irow, column=1, sticky=W, columnspan=2)
        self.browse = Button(Iframe, text="Browse", command=self.showfilebrowser)
        self.browse.grid(row=Irow, column=3, sticky=E)

        if Idescrip != "nohelp":
            self.help = Button(Iframe, text="?", command=self.showhelp)
            self.help.grid(row=Irow, column=4, sticky=E)



class FileSelectorSave:
    def enable(self):
        self.textentry.configure( state='normal')

    def disable(self):
        self.textentry.configure( state='disabled')

    def showhelp(self):
        helpWindows = SimpleDialog(self.root,
        text=self.descrip, buttons=["OK"], default=0, title="Help").go()

    def showfilebrowser(self):
        import tkFileDialog
        filetypes = [('Output File','*.out'), ('Any File','*.*')]
        file = tkFileDialog.asksaveasfilename(title='Save File',
                                              filetypes=filetypes)
        try:
            #Only do something if a file has been selected.
            if file:
                self.textentry.delete(0,255)
                self.textentry.insert(0, file)
        except:
            #If no file has been selected do nothing.
            pass

    def __init__(self, Iframe, Iroot, Irow,
        Iwidth, Istate, Ilabel, IstartText, Idescrip):
        self.text = StringVar()
        self.descrip = Idescrip
        self.root = Iroot
        self.textentry = Entry(Iframe, width=Iwidth, state=Istate, textvariable=self.text)
        self.textentry.insert(0, IstartText)
        self.label = Label(Iframe, text=Ilabel)
        self.label.grid(row=Irow, column=0, sticky=W)
        self.textentry.grid(row=Irow, column=1, sticky=W, columnspan=2)
        self.browse = Button(Iframe, text="Save as", command=self.showfilebrowser)
        self.browse.grid(row=Irow, column=3, sticky=E)

        if Idescrip != "nohelp":
            self.help = Button(Iframe, text="?", command=self.showhelp)
            self.help.grid(row=Irow, column=4, sticky=E)
