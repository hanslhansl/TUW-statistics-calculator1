import Pmw as PMW
import sys
import inspect
import tkinter as tk
from tkinter import ttk
import tkinter.font as font
import Funktionen
import tkintertable as tkt
import time
import numpy
import types



#Klasse für das Fenster
class main(tk.Tk):
    def __init__(self):
        self.distributions = {"Binomial" : None, "Normal" : None}
        self.categories = {"Konfidenzintervall" : None, "Hypothesentest" : None, "ANOVA" : None, "Verteilungen" : self.distributions, "Zwei Stichproben" : None}

        super().__init__()
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.mainFrame = tk.Frame(master=self)
        self.mainFrame.grid(sticky="NSWE")
        self.mainFrame.grid_rowconfigure(0, weight=1)
        self.mainFrame.grid_columnconfigure(1, weight=1)
        self.defaultFont = font.Font(family="Calibri", size=11)
        self.defaultFont2 = font.Font(family="Calibri", size=8)

        foo = inspect.getmembers(Funktionen, lambda member: inspect.isclass(member))
        self.functionsList = [x[1] for x in foo if x[1] != Funktionen.function and hasattr(x[1], "name")]
        self.functionsList = sorted(self.functionsList, key=lambda x: x.name)  #Sortiert alphabetisch
        self.functionsNamesList = [x.name for x in self.functionsList]

        self.l_dropdown = tk.Label(master=self.mainFrame, text="Available functions: ", font=self.defaultFont)
        self.l_dropdown.grid(row=0, column=0, sticky="NS")

        self.variableOptionMenu = tk.StringVar(self.mainFrame)
        self.variableOptionMenu.trace("w", self.chooseFunction)

        self.dropdown = tk.OptionMenu(self.mainFrame, self.variableOptionMenu, *self.functionsNamesList)
        self.dropdown.config(font=self.defaultFont)
        menu = self.nametowidget(self.dropdown.menuname)
        menu.config(font=self.defaultFont)
        self.dropdown.grid(row=0, column=1, sticky="WE")

        self.sep1 = ttk.Separator(master=self.mainFrame, orient="horizontal")
        self.sep1.grid(row=1, column=0, columnspan=2, sticky="EW")

        self.createMenuBar()

    #Erstellt die Menüleiste
    def createMenuBar(self):
        self.mainMenu = tk.Menu(self, tearoff=0)

        for name, content in self.categories.items():
            self._MenuBarHelper(parentMenu=self.mainMenu, key=name, value=content)

        self.config(menu=self.mainMenu)

    #Fügt der Menüleiste ein neues Untermenü hinzu
    def _MenuBarHelper(self, parentMenu, key, value):
        if value == None:
            menu = tk.Menu(parentMenu, tearoff=0)
            allClasses = [name for name in self.functionsNamesList if key in name]
            for name in allClasses:
                menu.add(itemType="command", label=name, command=lambda name=name: self.MenuBarCallback(entryName=name))

        elif type(value) == dict:
            menu = tk.Menu(parentMenu, tearoff=0)
            for name, content in value.items():
                self._MenuBarHelper(parentMenu=menu, key=name, value=content)

        else:
            raise Exception()

        parentMenu.add(itemType="cascade", label=key, menu=menu)


    def chooseFunction(self, *args):
        try:
            self.functionsFrame.destroy()
            del self.functionsFrame
        except AttributeError:
            pass
        self.functionsFrame = tk.Frame(master=self.mainFrame)
        self.functionsFrame.grid(row=2, column=0, columnspan=2, sticky="NSWE")
        self.functionsFrame.grid_columnconfigure(0, weight=1)
        self.functionsFrame.grid_columnconfigure(1, weight=1)

        index = self.functionsNamesList.index(self.variableOptionMenu.get())
        self.activeFunction = self.functionsList[index]()
        self.toolTip = PMW.Balloon(self.functionsFrame, initwait=200)

        #Beschreibung
        self.activeFunctionDescription = tk.Frame(master=self.functionsFrame)
        text = self.activeFunction.description
        if text != "":
            foo = tk.Text(self.activeFunctionDescription, width=200, height=8, borderwidth=0, font=self.defaultFont, relief=tk.FLAT)#, background=self.resultFrame.cget("background"))
            foo.tag_configure("sub", offset=-3, font=self.defaultFont2)
            scrollb = ttk.Scrollbar(self.activeFunctionDescription, command=foo.yview)
            while True:
                a = text.find("<sub>")
                b = text.find("<\sub>")
                if a != -1 and b != -1:
                    first, sec = text.split("<sub>", 1)
                    sec, text = sec.split("<\sub>", 1)
                    foo.insert("insert", first, "", sec, "sub")
                elif a == -1 and b == -1:
                    foo.insert("insert", text)
                    break
                else:
                    raise Exception(f"Unbalanced tag placement: {text}")
            foo.configure(state="disabled")
            text = ""

            foo.grid(row=0, column=0)
            scrollb.grid(row=0, column=1, sticky='ns')
            foo['yscrollcommand'] = scrollb.set
        self.activeFunctionDescription.grid(row=0, column=0, columnspan=2)

        self.radioFrame = tk.Frame(master=self.functionsFrame)
        self.radioFrame.grid(row=1, column=0, columnspan=2, sticky="NSWE")
        self.radioLabels = []
        self.radioButtons = []
        tempDict = self.activeFunction.necessaryValues
        self.radioVariable = tk.StringVar()
        i = 0
        for key in tempDict:
            self.radioButtons.append(tk.Radiobutton(master=self.radioFrame, text=key, variable=self.radioVariable, value=key, command=self.radioCallback, tristatevalue=0))
            self.radioButtons[-1].grid(row=0, column=i)
            self.radioFrame.grid_columnconfigure(i, weight=1)
            i += 1

        self.sep2 = ttk.Separator(master=self.functionsFrame, orient="horizontal")
        self.sep2.grid(row=2, column=0, columnspan=2, sticky="EW")

        self.inputFrame = tk.Frame(master=self.functionsFrame)
        self.inputFrame.grid(row=3, column=0, columnspan=2, sticky="NSWE")
        self.inputDict = {}

        self.resultFrame = tk.Frame(master=self.functionsFrame)
        self.resultFrame.grid(row=4, column=0, columnspan=2, sticky="NSWE")

        if type(self.activeFunction.variableNames) == list:
            pass
        else:
            raise Exception()

        i = 0
        for key in self.activeFunction.variableNames:
            value = self.activeFunction.extendedExplanations[key]
            if type(value) == list:
                frame = tk.Frame(master=self.inputFrame)
                frame.grid(row=0, column=2*i)
                self.inputFrame.grid_columnconfigure(2*i, weight=50)
                foo = {x : "" for x in value}
                data = {"1" : foo}
                cols = len(value)
                table = tkt.TableCanvas(frame, data=data, width=150*cols, rowheaderwidth=20, cols=cols, cellwidth=150, height=90)
                table.injection = types.MethodType(smallInjection, table)
                table.show()

                self.inputDict[key] = [table.model]

            elif type(key) == str:
                label = tk.Label(master=self.inputFrame, text=key+": ", font=self.defaultFont)
                label.grid(row=0, column=2*i)
                self.toolTip.bind(label, str(value))

                var = tk.StringVar(master=self.inputFrame)
                entry = tk.Entry(master=self.inputFrame, textvariable=var, font=self.defaultFont, state=tk.DISABLED, width=15)
                entry.grid(row=0, column=2*i+1)
                self.inputFrame.grid_columnconfigure(2*i+1, weight=1)

                self.inputDict[key] = [label, entry, var]

            i += 1

        self.b_calculate = tk.Button(master=self.functionsFrame, text="Calculate", command=self.calculate, font=self.defaultFont, state=tk.DISABLED)
        self.b_calculate.grid(row=5, column=0)

        self.b_clear = tk.Button(master=self.functionsFrame, text="Clear", command=self.chooseFunction, font=self.defaultFont)
        self.b_clear.grid(row=5, column=1)


    def calculate(self):
        tempDict = {}
        for key, value in self.inputDict.items():
            if len(value) == 3:
                tempDict[key] = value[2].get()
            elif len(value) == 1:
                tempDict[key] = value[0]

        try:
            self.resultFrame.destroy()
            del self.resultFrame
        except AttributeError:
            pass
        self.result = self.activeFunction(tempDict)
        self.resultFrame = tk.Frame(master=self.functionsFrame)
        self.resultFrame.grid(row=4, column=0, columnspan=2, sticky="NSWE")

        text = ""
        i = 0
        for key, value in self.result.items():
            if type(value) in (dict, tuple, str, int, float, numpy.float64) or key == list(self.result.keys())[-1]:

                if type(value) in (str, int, float, numpy.float64, tuple):
                    if type(value) == tuple:
                        if len(value) != 1:
                            raise Exception()
                        text = text[:-3]
                        text = text + f"\n{key}:\n{self.result[key][0]}\n"
                    else:
                        text = text + f"{key} = {self.result[key]};  "

                if type(value) in (dict, ) or key == list(self.result.keys())[-1]:
                    text = text[:-3]

                    if text != "":
                        foo = tk.Text(self.resultFrame, width=200, height=10, borderwidth=0, font=self.defaultFont, relief=tk.FLAT)#, background=self.resultFrame.cget("background"))
                        foo.tag_configure("sub", offset=-3, font=self.defaultFont2)
                        scrollb = ttk.Scrollbar(self.resultFrame, command=foo.yview)

                        while True:
                            a = text.find("<sub>")
                            b = text.find("<\sub>")
                            if a != -1 and b != -1:
                                first, sec = text.split("<sub>", 1)
                                sec, text = sec.split("<\sub>", 1)
                                foo.insert("insert", first, "", sec, "sub")
                            elif a == -1 and b == -1:
                                foo.insert("insert", text)
                                break
                            else:
                                raise Exception(f"Unbalanced tag placement: {text}")
                        foo.configure(state="disabled")
                        text = ""

                        foo.grid(row=i, column=0)
                        scrollb.grid(row=i, column=1, sticky='ns')
                        foo['yscrollcommand'] = scrollb.set
                        i+=1

                    if type(value) == dict:
                        frame = tk.Frame(master=self.resultFrame)
                        self.resultFrame.columnconfigure(0, weight=1)
                        frame.grid(row=i, column=0, sticky="NSWE")
                        rows = len(value)
                        table = tkt.TableCanvas(frame, data=value, rows=rows, height=20*rows)
                        table.adjustColumnWidths()
                        table.show()
                        i+=1
            else:
                raise Exception(f"Type of value is {type(value)}")


    def _report_callback_exception(self, *args, **kwaargs):
        print("============================")
        print("Error Handler:")
        print("args: ", args)
        print("kwaargs: ", kwaargs)
        print("============================")
        raise args[1]


    def radioCallback(self, *args):
        not_disable = self.activeFunction.necessaryValues[self.radioVariable.get()]

        self.inputFrame.destroy()
        del self.inputFrame
        self.inputFrame = tk.Frame(master=self.functionsFrame)
        self.inputFrame.grid(row=3, column=0, columnspan=2, sticky="NSWE")
        self.inputDict = {}

        self.resultFrame.destroy()
        del self.resultFrame
        self.resultFrame = tk.Frame(master=self.functionsFrame)
        self.resultFrame.grid(row=4, column=0, columnspan=2, sticky="NSWE")

        i = 0
        for key in not_disable:
            value = self.activeFunction.extendedExplanations[key]
            if type(value) == list:
                if type(value[0]) == str:
                    frame = tk.Frame(master=self.inputFrame)
                    frame.grid(row=0, column=2*i, sticky="WE")
                    self.inputFrame.grid_columnconfigure(2*i, weight=50)
                    foo = {x : "" for x in value}
                    data = {"1" : foo}
                    cols = len(value)
                    table = tkt.TableCanvas(frame, data=data, width=150*cols, height=200, cellwidth=150, rowheaderwidth=15, cols=cols, rows=1)
                    table.injection = types.MethodType(smallInjection, table)

                if type(value[0]) == list:
                    value = value[0]
                    if len(value) != 2:
                        raise Exception()

                    frame = tk.Frame(master=self.inputFrame)
                    frame.grid(row=0, column=2*i, sticky="WE")
                    self.inputFrame.grid_columnconfigure(2*i, weight=50)
                    foo = {f"1. {value[0]}" : ""}
                    data = {f"{value[1]}" : foo}
                    table = tkt.TableCanvas(frame, data=data, width=150, height=200, cellwidth=150, rowheaderwidth=15, cols=1, rows=1)
                    table.injection = types.MethodType(bigInjection, table)

                table.show()
                self.inputDict[key] = [table.model]

            elif type(key) == str:
                label = tk.Label(master=self.inputFrame, text=key+": ", font=self.defaultFont)
                label.grid(row=0, column=2*i, sticky="WE")
                self.toolTip.bind(label, str(value))

                var = tk.StringVar(master=self.inputFrame)
                entry = tk.Entry(master=self.inputFrame, textvariable=var, font=self.defaultFont, width=15)
                entry.grid(row=0, column=2*i+1, sticky="WE")
                self.inputFrame.grid_columnconfigure(2*i+1, weight=1)

                self.inputDict[key] = [label, entry, var]

            i += 1

        self.b_calculate["state"] = tk.NORMAL


    def MenuBarCallback(self, entryName):
        self.variableOptionMenu.set(entryName)


def fakeClick(self, row=None, col=None):
    if row == None:
        row = self.getSelectedRow()
    if col == None:
        col = self.getSelectedColumn()

    #which row and column is the click inside?
    self.clearSelected()
    self.allrows = False
    rowclicked = row
    colclicked = col
    self.focus_set()
    if self.mode == 'formula':
        self.handleFormulaClick(rowclicked, colclicked)
        return
    if hasattr(self, 'cellentry'):
        self.cellentry.destroy()
    #ensure popup menus are removed if present
    if hasattr(self, 'rightmenu'):
        self.rightmenu.destroy()
    if hasattr(self.tablecolheader, 'rightmenu'):
        self.tablecolheader.rightmenu.destroy()

    self.startrow = rowclicked
    self.endrow = rowclicked
    self.startcol = colclicked
    self.endcol = colclicked
    #reset multiple selection list
    self.multiplerowlist=[]
    self.multiplerowlist.append(rowclicked)
    if rowclicked is None or colclicked is None:
        return
    if self.read_only is True:
        return
    if 0 <= rowclicked < self.rows and 0 <= colclicked < self.cols:
        self.setSelectedRow(rowclicked)
        self.setSelectedCol(colclicked)
        self.drawSelectedRect(self.currentrow, self.currentcol)
        self.drawSelectedRow()
        self.tablerowheader.drawSelectedRows(rowclicked)
        coltype = self.model.getColumnType(colclicked)
        if coltype == 'text' or coltype == 'number':
            self.drawCellEntry(rowclicked, colclicked)
    return


def handle_left_click(self, event):
    """Respond to a single press"""

    #which row and column is the click inside?
    self.clearSelected()
    self.allrows = False
    rowclicked = self.get_row_clicked(event)
    colclicked = self.get_col_clicked(event)

    self.focus_set()
    if self.mode == 'formula':
        self.handleFormulaClick(rowclicked, colclicked)
        return
    if hasattr(self, 'cellentry'):
        self.cellentry.destroy()
    #ensure popup menus are removed if present
    if hasattr(self, 'rightmenu'):
        self.rightmenu.destroy()
    if hasattr(self.tablecolheader, 'rightmenu'):
        self.tablecolheader.rightmenu.destroy()

    self.startrow = rowclicked
    self.endrow = rowclicked
    self.startcol = colclicked
    self.endcol = colclicked
    #reset multiple selection list
    self.multiplerowlist=[]
    self.multiplerowlist.append(rowclicked)
    if rowclicked is None or colclicked is None:
        return
    if self.read_only is True:
        return
    if 0 <= rowclicked < self.rows and 0 <= colclicked < self.cols:
        self.setSelectedRow(rowclicked)
        self.setSelectedCol(colclicked)
        self.drawSelectedRect(self.currentrow, self.currentcol)
        self.drawSelectedRow()
        self.tablerowheader.drawSelectedRows(rowclicked)
        coltype = self.model.getColumnType(colclicked)
        if coltype == 'text' or coltype == 'number':
            self.drawCellEntry(rowclicked, colclicked)

    self.injection(rowclicked, colclicked)
    return


def smallInjection(self, rowclicked, colclicked):
    if rowclicked == self.rows-1:
        if not all([self.model.getCellRecord(self.model.getRowCount()-1, k) in (None, "") for k in range(self.model.getColumnCount())]):
            self.addRow()
            self.fakeClick(row=self.rows-2)
    elif rowclicked < self.rows-1:
        if all([self.model.getCellRecord(self.model.getRowCount()-1, k) in (None, "") for k in range(self.model.getColumnCount())]):
            self.model.deleteRow(self.model.getRowCount()-1)
            self.fakeClick()
            self.redrawTable()


def bigInjection(self, rowclicked, colclicked):
    if rowclicked >= self.rows-1 or colclicked >= self.cols-1:
        if rowclicked <= self.rows-1 and colclicked <= self.cols-1:
            if rowclicked == self.rows-1:
                if not all([self.model.getCellRecord(self.model.getRowCount()-1, k) in (None, "") for k in range(self.model.getColumnCount())]):
                    self.addRow()
                    self.fakeClick(row=self.rows-2)

            if colclicked == self.cols-1:
                if not all([self.model.getCellRecord(k, self.model.getColumnCount()-1) in (None, "") for k in range(self.model.getRowCount())]):
                    lastcolname = self.model.getColumnLabel(self.cols-1)
                    no, name = lastcolname.split(". ", 1)
                    no = int(no)
                    self.addColumn(f"{no+1}. {name}")
                    self.fakeClick(col=self.cols-2)

    elif colclicked < self.cols-1 and rowclicked < self.rows-1:
        allRow = all([self.model.getCellRecord(self.model.getRowCount()-1, k) in (None, "") for k in range(self.model.getColumnCount())])
        allCol = all([self.model.getCellRecord(k, self.model.getColumnCount()-1) in (None, "") for k in range(self.model.getRowCount())])
        if allRow:
            self.model.deleteRow(self.model.getRowCount()-1)
        if allCol:
            self.model.deleteColumn(self.model.getColumnCount()-1)
        self.fakeClick()
        self.redrawTable()


tkt.TableCanvas.drawTooltip = lambda *args, **kwargs: None
tkt.TableCanvas.handle_left_click = handle_left_click
tkt.TableCanvas.fakeClick = fakeClick
tkt.TableCanvas.injection = lambda *args, **kwargs: None

root = main()
#root.geometry("1500x700")
root.minsize(800, 100)
root.title("Stochastik Rechner")


root.mainloop()
