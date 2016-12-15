from PyQt4 import QtCore, QtGui, QtWebKit, QtNetwork
import qdarkstyle
import pyqtgraph as pg
from obspy.clients.fdsn.client import Client
from obspy.core import UTCDateTime
from obspy.geodetics import degrees2kilometers
import sys
import pandas as pd
import numpy as np
from DateAxisItem import DateAxisItem
import math


class PandasModel(QtCore.QAbstractTableModel):
    """
    Class to populate a table view with a pandas dataframe
    """

    def __init__(self, data, isc_nm=None, oth_nm=None, parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self._data = np.array(data.values)
        self._cols = data.columns
        self.r, self.c = np.shape(self._data)

        # Column headers for tables
        self.isc_col_header = ['isc_ind', 'Event ID', 'qtime', 'Lat (dd)', 'Lon  (dd)', 'Depth (km)', 'Mag']
        self.oth_col_header = ['isc_ind', 'Event ID', 'qtime', 'Lat (dd)', 'Lon  (dd)', 'Depth (km)', 'Mag']
        self.matched_col_header = ['oth_ind', 'Event ID', 'qtime', 'Lat (dd)', 'Lon (dd)', 'Depth (km)', 'Mag',
                                   'isc_ind', 'Event ID Match', 'qtime_match', 'Lat Match (dd)', 'Lon Match (dd)',
                                   'Depth Match (km)', 'Mag Match', 'Time Diff (s)', 'Distance Diff (km)',
                                   'Depth Diff (km)', 'Mag Diff']

        self.isc_nm = isc_nm
        self.oth_nm = oth_nm

        self.matched_color = QtGui.QColor(0, 0, 0)
        self.matched_color.setNamedColor('#008000')
        self.matched_color.setAlpha(100)
        self.isc_nm_color = QtGui.QColor(0, 0, 0)
        self.isc_nm_color.setNamedColor('#4169E1')
        self.isc_nm_color.setAlpha(100)
        self.oth_nm_color = QtGui.QColor(0, 0, 0)
        self.oth_nm_color.setNamedColor('#FFA500')
        self.oth_nm_color.setAlpha(100)

    def rowCount(self, parent=None):
        return self.r

    def columnCount(self, parent=None):
        return self.c

    def data(self, index, role=QtCore.Qt.DisplayRole):

        if index.isValid():
            if role == QtCore.Qt.DisplayRole:
                return self._data[index.row(), index.column()]
            if role == QtCore.Qt.BackgroundRole:
                if not self.isc_nm == None:
                    bool_in = index.row() in self.isc_nm
                    if bool_in:
                        return QtCore.QVariant(QtGui.QBrush(self.isc_nm_color))
                    elif not bool_in:
                        return QtCore.QVariant(QtGui.QBrush(self.matched_color))
                if not self.oth_nm == None:
                    bool_in = index.row() in self.oth_nm
                    if bool_in:
                        return QtCore.QVariant(QtGui.QBrush(self.oth_nm_color))
                    elif not bool_in:
                        return QtCore.QVariant(QtGui.QBrush(self.matched_color))

                if self.isc_nm == None and self.oth_nm == None:
                    return QtCore.QVariant(QtGui.QBrush(self.matched_color))
        return None

    def headerData(self, p_int, orientation, role):
        if role == QtCore.Qt.DisplayRole:
            if orientation == QtCore.Qt.Horizontal:
                if not self.isc_nm == None:
                    return self.isc_col_header[p_int]
                elif not self.oth_nm == None:
                    return self.oth_col_header[p_int]
                elif self.isc_nm == None and self.oth_nm == None:
                    return self.matched_col_header[p_int]
            elif orientation == QtCore.Qt.Vertical:
                return p_int
        return None


class TableDialog(QtGui.QDialog):
    """
    Class to create a separate child window to display the isc, local and matched Pandas arrays as 3 tables
    """

    def __init__(self, parent=None, isc_df=None, oth_df=None, matched_df=None, isc_nm=None, oth_nm=None):
        super(TableDialog, self).__init__(parent)

        self.isc_df = isc_df
        self.oth_df = oth_df
        self.matched_df = matched_df

        self.isc_nm = isc_nm
        self.oth_nm = oth_nm

        self.initUI()

    def initUI(self):
        self.layout = QtGui.QVBoxLayout(self)

        self.isc_event_table_view = QtGui.QTableView()
        self.oth_event_table_view = QtGui.QTableView()
        self.matched_event_table_view = QtGui.QTableView()

        self.isc_event_table_view.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.oth_event_table_view.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.matched_event_table_view.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)

        self.layout.addWidget(self.isc_event_table_view)
        self.layout.addWidget(self.oth_event_table_view)
        self.layout.addWidget(self.matched_event_table_view)

        self.setLayout(self.layout)

        self.isc_model = PandasModel(self.isc_df, isc_nm=self.isc_nm)
        self.oth_model = PandasModel(self.oth_df, oth_nm=self.oth_nm)
        self.matched_model = PandasModel(self.matched_df)

        self.isc_event_table_view.setModel(self.isc_model)
        self.oth_event_table_view.setModel(self.oth_model)
        self.matched_event_table_view.setModel(self.matched_model)

        self.setWindowTitle('Tables')
        self.show()


class MainWindow(QtGui.QWidget):
    """
    Main Window for compare_catalogues GUI
    """

    def __init__(self):
        super(MainWindow, self).__init__()
        self.setupUi()
        self.show()
        self.raise_()
        QtGui.QApplication.instance().focusChanged.connect(self.changed_widget_focus)

        self.get_catalogues()

        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)

    def setupUi(self):

        vbox = QtGui.QVBoxLayout()
        self.setLayout(vbox)

        hbox = QtGui.QHBoxLayout()

        self.compare_catalogues_button = QtGui.QPushButton('Compare Catalogues!')
        self.compare_catalogues_button.released.connect(self.compare_catalogues)
        hbox.addWidget(self.compare_catalogues_button)

        self.icon_size_by_depth_button = QtGui.QPushButton('Icon Size by Depth')
        self.icon_size_by_depth_button.released.connect(self.icon_size_by_depth)
        self.icon_size_by_depth_button.setEnabled(False)
        hbox.addWidget(self.icon_size_by_depth_button)

        self.icon_size_by_mag_button = QtGui.QPushButton('Icon Size by Mag')
        self.icon_size_by_mag_button.released.connect(self.icon_size_by_mag)
        self.icon_size_by_mag_button.setEnabled(False)
        hbox.addWidget(self.icon_size_by_mag_button)

        self.icon_size_reset_button = QtGui.QPushButton('Icon Size Reset')
        self.icon_size_reset_button.released.connect(self.icon_size_reset)
        self.icon_size_reset_button.setEnabled(False)
        hbox.addWidget(self.icon_size_reset_button)

        vbox.addLayout(hbox)

        web_view = self.web_view = QtWebKit.QWebView()

        web_view.page().mainFrame().addToJavaScriptWindowObject("MainWindow", self)
        web_view.page().setLinkDelegationPolicy(QtWebKit.QWebPage.DelegateAllLinks)
        web_view.load(QtCore.QUrl('map.html'))
        web_view.loadFinished.connect(self.onLoadFinished)
        web_view.linkClicked.connect(QtGui.QDesktopServices.openUrl)

        vbox.addWidget(web_view)

        self.graph_view = pg.GraphicsLayoutWidget()

        vbox.addWidget(self.graph_view)

    def onLoadFinished(self):
        with open('map.js', 'r') as f:
            frame = self.web_view.page().mainFrame()
            frame.evaluateJavaScript(f.read())

    def changed_widget_focus(self):
        # deselect the highlighted point in the scatter graphs
        # if this is not done then the event will be continuously highlighted in the map view
        # even if another event is clicked on via the map view
        try:
            if not QtGui.QApplication.focusWidget() == self.graph_view:
                self.scatter_point_deselect()
        except AttributeError:
            pass

    def get_catalogues(self):
        isc_catalogue = Client("http://isc-mirror.iris.washington.edu")
        oth_catalogue = Client('http://sismo.iag.usp.br')

        t_start = UTCDateTime('2016-12-01T00:00:00')
        t_end = UTCDateTime('2016-12-02T00:00:00')

        print('\nRequesting Earthquake Catalogues from Remote Servers.....')

        self.isc_cat = isc_catalogue.get_events(starttime=t_start, endtime=t_end)
        self.oth_cat = oth_catalogue.get_events(starttime=t_start, endtime=t_end)

        # create empty data frame
        self.isc_df = pd.DataFrame(data=None, columns=['isc_ind', 'event_id', 'qtime', 'lat', 'lon', 'depth', 'mag'])
        self.oth_df = pd.DataFrame(data=None, columns=['oth_ind', 'event_id', 'qtime', 'lat', 'lon', 'depth', 'mag'])

        # iterate through the events in oth cat
        for _i, event in enumerate(self.oth_cat):
            print "\r     Parsing event from Local Cat", _i + 1, ' of ', len(self.oth_cat), ' ....',
            sys.stdout.flush()

            # Get quake origin info
            origin_info = event.preferred_origin() or event.origins[0]

            try:
                mag_info = event.preferred_magnitude() or event.magnitudes[0]
                magnitude = mag_info.mag
            except IndexError:
                # No magnitude for event
                magnitude = None

            self.oth_df.loc[_i] = [int(_i), str(event.resource_id.id), int(origin_info.time.timestamp),
                                   origin_info.latitude, origin_info.longitude,
                                   origin_info.depth, magnitude]

        # Convert the depth field to Km from m
        self.oth_df.loc[:, 'depth'] *= (1.0 / 1000)

        print('\n')
        # iterate through the events in isc cat
        for _i, event in enumerate(self.isc_cat):
            print "\r     Parsing event from ISC Cat", _i + 1, ' of ', len(self.isc_cat), ' ....',
            sys.stdout.flush()
            # Get quake origin info
            origin_info = event.preferred_origin() or event.origins[0]
            try:
                mag_info = event.preferred_magnitude() or event.magnitudes[0]
                magnitude = mag_info.mag
            except IndexError:
                # No magnitude for event
                magnitude = None
            self.isc_df.loc[_i] = [int(_i), str(event.resource_id.id), int(origin_info.time.timestamp),
                                   origin_info.latitude, origin_info.longitude,
                                   origin_info.depth, magnitude]

        # Convert the depth field to Km from m
        self.isc_df.loc[:, 'depth'] *= (1.0 / 1000)

        # =====================Finding matching events =======================

        print('\nFinding Matching Events.....')

        global match_index
        global length_oth_df
        match_index = 0
        length_oth_df = len(self.oth_df)

        def get_isc_match(row):
            global match_index
            global length_oth_df
            print "\r     Matching event from Local Cat", match_index, ' of ', length_oth_df, ' ....',
            sys.stdout.flush()
            temp = self.isc_df_drop.apply(lambda x: abs(x - row), axis=1)  # Pandas DF
            # NaNs are treated as small
            smallest_temp = temp.nsmallest(2, columns=['lat', 'lon', 'qtime']).iloc[0]  # Pandas Series

            distance_diff = degrees2kilometers(math.sqrt(abs(smallest_temp['lat']) ** 2 +
                                                         abs(smallest_temp['lon']) ** 2))

            isc_index = smallest_temp.name

            if smallest_temp['qtime'] <= 15 and \
                    (abs(smallest_temp['lon']) <= 1 or np.isnan(smallest_temp['lon'])) and \
                    (abs(smallest_temp['lat']) <= 1 or np.isnan(smallest_temp['lat'])):
                ret_s = pd.Series([isc_index, self.isc_df.loc[isc_index, 'event_id'],
                                   self.isc_df.loc[isc_index, 'qtime'], self.isc_df.loc[isc_index, 'lat'],
                                   self.isc_df.loc[isc_index, 'lon'],
                                   self.isc_df.loc[isc_index, 'depth'], self.isc_df.loc[isc_index, 'mag'],
                                   smallest_temp['qtime'], distance_diff, smallest_temp['depth'], smallest_temp['mag']],
                                  index=['isc_ind', 'event_id_match', 'qtime_match', 'lat_match', 'lon_match',
                                         'depth_match', 'mag_match', 'qtime_diff', 'dist_diff', 'depth_diff',
                                         'mag_diff'])
            else:
                ret_s = pd.Series([None, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  index=['isc_ind', 'event_id_match', 'qtime_match', 'lat_match', 'lon_match',
                                         'depth_match', 'mag_match', 'qtime_diff', 'dist_diff', 'depth_diff',
                                         'mag_diff'])

            match_index += 1
            return ret_s

        # Drop the event_id column (strings) from the data frame to apply vectorised function
        self.oth_df_drop = self.oth_df.drop('event_id', axis=1)
        self.isc_df_drop = self.isc_df.drop('event_id', axis=1)
        self.matched_df = pd.concat((self.oth_df, self.oth_df_drop.apply(get_isc_match, axis=1)), axis=1)

        # drop a row from the matched df if isc_ind in matched df is NaN
        # (I.e. there was no matching earthquake in isc cat)
        self.matched_df.dropna(subset=['isc_ind'], inplace=True)
        self.matched_df.reset_index(drop=True, inplace=True)

        # find isc events not matched and oth events not matched
        self.isc_not_matched_df = self.isc_df[~self.isc_df['isc_ind'].isin(self.matched_df['isc_ind'])]
        self.oth_not_matched_df = self.oth_df[~self.oth_df['oth_ind'].isin(self.matched_df['oth_ind'])]

    def compare_catalogues(self):
        self.compare_catalogues_button.setEnabled(False)
        self.icon_size_by_depth_button.setEnabled(True)
        self.icon_size_by_mag_button.setEnabled(True)
        self.icon_size_reset_button.setEnabled(True)

        # Plot the matched events
        for row_index, row in self.matched_df.iterrows():
            js_call = "addEvent('{event_id}', '{df_id}', {row_index}, " \
                      "{latitude}, {longitude}, '{a_color}', '{p_color}');" \
                .format(event_id=row['event_id'], df_id="matched", row_index=int(row_index), latitude=row['lat'],
                        longitude=row['lon'], a_color="Red",
                        p_color="#008000")
            self.web_view.page().mainFrame().evaluateJavaScript(js_call)

        # Plot ISC events (not matched)
        for row_index, row in self.isc_not_matched_df.iterrows():
            js_call = "addEvent('{event_id}', '{df_id}', {row_index}, " \
                      "{latitude}, {longitude}, '{a_color}', '{p_color}');" \
                .format(event_id=row['event_id'], df_id="isc", row_index=int(row_index), latitude=row['lat'],
                        longitude=row['lon'], a_color="Red",
                        p_color="#4169E1")
            self.web_view.page().mainFrame().evaluateJavaScript(js_call)

        # Plot oth events (not matched)
        for row_index, row in self.oth_not_matched_df.iterrows():
            js_call = "addEvent('{event_id}', '{df_id}', {row_index}, " \
                      "{latitude}, {longitude}, '{a_color}', '{p_color}');" \
                .format(event_id=row['event_id'], df_id="oth", row_index=int(row_index), latitude=row['lat'],
                        longitude=row['lon'], a_color="Red",
                        p_color="#FFA500")
            self.web_view.page().mainFrame().evaluateJavaScript(js_call)

        js_call = "createOverlay();"
        self.web_view.page().mainFrame().evaluateJavaScript(js_call)

        self.table_accessor = None

        isc_nm, oth_nm = self.table_view_create_color()

        self.tbld = TableDialog(parent=self, isc_df=self.isc_df, oth_df=self.oth_df,
                                matched_df=self.matched_df, isc_nm=isc_nm, oth_nm=oth_nm)

        # Lookup Dictionary for table views
        self.tbl_view_dict = {"isc": self.tbld.isc_event_table_view,
                              "oth": self.tbld.oth_event_table_view,
                              "matched": self.tbld.matched_event_table_view}

        # Create a new table_accessor dictionary for this class
        self.table_accessor = {self.tbld.isc_event_table_view: [self.isc_df, range(0, len(self.isc_df))],
                               self.tbld.oth_event_table_view: [self.oth_df, range(0, len(self.oth_df))],
                               self.tbld.matched_event_table_view: [self.matched_df, range(0, len(self.matched_df))]}

        self.tbld.isc_event_table_view.clicked.connect(self.table_view_clicked)
        self.tbld.oth_event_table_view.clicked.connect(self.table_view_clicked)
        self.tbld.matched_event_table_view.clicked.connect(self.table_view_clicked)

        self.tbld.isc_event_table_view.horizontalHeader().sectionClicked.connect(self.table_view_headerclicked)
        self.tbld.oth_event_table_view.horizontalHeader().sectionClicked.connect(self.table_view_headerclicked)
        self.tbld.matched_event_table_view.horizontalHeader().sectionClicked.connect(self.table_view_headerclicked)

        # Double click and zoom is making markers disappear when clicking on them
        # self.tbld.isc_event_table_view.doubleClicked.connect(self.table_view_dblclicked)
        # self.tbld.oth_event_table_view.doubleClicked.connect(self.table_view_dblclicked)
        # self.tbld.matched_event_table_view.doubleClicked.connect(self.table_view_dblclicked)

        self.scatter_build_graph()

    def icon_size_by_depth(self):
        # Change the sizes of the icons in map view to represent depth difference between matched events
        depth_norm_series = self.matched_df.loc[:, 'depth_diff'] \
            .apply(lambda x: (x / self.matched_df.loc[:, 'depth_diff'].max()))
        depth_size_list = (depth_norm_series + 0.3).multiply(18, fill_value=0.1).tolist()

        js_call = "arrange_marker_by_depth({sizes});" \
            .format(sizes=depth_size_list)
        self.web_view.page().mainFrame().evaluateJavaScript(js_call)

    def icon_size_by_mag(self):
        # Change the sizes of the icons in map view to represent magnitude difference between matched events
        mag_norm_series = self.matched_df.loc[:, 'mag_diff'] \
            .apply(lambda x: (x / self.matched_df.loc[:, 'mag_diff'].max()))
        mag_size_list = (mag_norm_series + 0.3).multiply(18, fill_value=0.1).tolist()

        js_call = "arrange_marker_by_depth({sizes});" \
            .format(sizes=mag_size_list)
        self.web_view.page().mainFrame().evaluateJavaScript(js_call)

    def icon_size_reset(self):
        js_call = "resetMarkerSize();"
        self.web_view.page().mainFrame().evaluateJavaScript(js_call)
        js_call = "setAllInactive();"
        self.web_view.page().mainFrame().evaluateJavaScript(js_call)

    def table_view_highlight(self, focus_widget, row_index):

        if focus_widget == self.tbld.matched_event_table_view:
            selected_row = self.matched_df.loc[row_index]

            # get the index of the the ISC and oth matched quake
            isc_index = selected_row['isc_ind']
            oth_index = selected_row['oth_ind']

            # Find the row_numbers of these two indexes
            matched_row_number = self.table_accessor[self.tbld.matched_event_table_view][1].index(row_index)
            isc_row_number = self.table_accessor[self.tbld.isc_event_table_view][1].index(isc_index)
            oth_row_number = self.table_accessor[self.tbld.oth_event_table_view][1].index(oth_index)

            # Highlight/Select the  rows in the table
            self.tbld.matched_event_table_view.selectRow(matched_row_number)
            self.tbld.isc_event_table_view.selectRow(isc_row_number)
            self.tbld.oth_event_table_view.selectRow(oth_row_number)

            # Highlight the marker on the map
            js_call = "highlightEvent('{event_id}');".format(event_id=selected_row['event_id'])
            self.web_view.page().mainFrame().evaluateJavaScript(js_call)

        elif focus_widget == self.tbld.isc_event_table_view:
            self.tbld.oth_event_table_view.clearSelection()
            self.tbld.matched_event_table_view.clearSelection()

            selected_row = self.isc_df.loc[row_index]

            # Find the row_number of this index
            isc_row_number = self.table_accessor[focus_widget][1].index(row_index)
            focus_widget.selectRow(isc_row_number)

            # Highlight the marker on the map
            js_call = "highlightEvent('{event_id}');".format(event_id=selected_row['event_id'])
            self.web_view.page().mainFrame().evaluateJavaScript(js_call)

        elif focus_widget == self.tbld.oth_event_table_view:
            self.tbld.isc_event_table_view.clearSelection()
            self.tbld.matched_event_table_view.clearSelection()

            selected_row = self.oth_df.loc[row_index]

            # Find the row_number of this index
            oth_row_number = self.table_accessor[focus_widget][1].index(row_index)
            focus_widget.selectRow(oth_row_number)

            # Highlight the marker on the map
            js_call = "highlightEvent('{event_id}');".format(event_id=selected_row['event_id'])
            self.web_view.page().mainFrame().evaluateJavaScript(js_call)

    def table_view_create_color(self):
        # Find the row number in the table_view corresponding to events that are matched,
        # not_matched (local cat), not matched (international cat)
        find = lambda searchList, elem: [[i for i, x in enumerate(searchList) if x == e][0] for e in elem]

        # The first time initializing the table
        if self.table_accessor == None:
            if not self.oth_not_matched_df.empty:
                # first unmatched events in local catalogue
                oth_not_matched_rn_list = self.oth_not_matched_df.loc[:, 'oth_ind'].tolist()
            else:
                oth_not_matched_rn_list = None
            if not self.isc_not_matched_df.empty:
                # now unmatched events in local catalogue
                isc_not_matched_rn_list = self.isc_not_matched_df.loc[:, 'isc_ind'].tolist()
            else:
                isc_not_matched_rn_list = None

        else:
            if not self.oth_not_matched_df.empty:
                # first unmatched events in local catalogue
                oth_not_matched_rn_list = find(self.table_accessor[self.tbld.oth_event_table_view][1],
                                               self.oth_not_matched_df.loc[:, 'oth_ind'].tolist())
            else:
                oth_not_matched_rn_list = None
            if not self.isc_not_matched_df.empty:
                # now unmatched events in local catalogue
                isc_not_matched_rn_list = find(self.table_accessor[self.tbld.isc_event_table_view][1],
                                               self.isc_not_matched_df.loc[:, 'isc_ind'].tolist())
            else:
                isc_not_matched_rn_list = None

        return (isc_not_matched_rn_list, oth_not_matched_rn_list)

    def table_view_headerclicked(self, logicalIndex):
        focus_widget = QtGui.QApplication.focusWidget()
        table_df = self.table_accessor[focus_widget][0]

        header = focus_widget.horizontalHeader()

        self.order = header.sortIndicatorOrder()
        table_df.sort_values(by=table_df.columns[logicalIndex],
                             ascending=self.order, inplace=True)

        self.table_accessor[focus_widget][1] = table_df.index.tolist()

        isc_nm, oth_nm = self.table_view_create_color()

        if focus_widget == self.tbld.isc_event_table_view:
            self.model = PandasModel(table_df, isc_nm=isc_nm)
        elif focus_widget == self.tbld.oth_event_table_view:
            self.model = PandasModel(table_df, oth_nm=oth_nm)
        else:
            self.model = PandasModel(table_df)

        focus_widget.setModel(self.model)
        focus_widget.update()

    def table_view_dblclicked(self):
        self.table_view_clicked()
        self.panMap(self.selected_row['lon'], self.selected_row['lat'])

    def table_view_clicked(self):
        focus_widget = QtGui.QApplication.focusWidget()
        row_number = focus_widget.selectionModel().selectedRows()[0].row()
        row_index = self.table_accessor[focus_widget][1][row_number]
        # Highlight/Select the current row in the table
        self.table_view_highlight(focus_widget, row_index)

    def scatter_build_graph(self):
        # Plot the Magnitude Difference of matched events, Depth Difference of matched events and the
        # Unmatched events as scatter plots
        self.graph_view.clear()

        # Set the colors of the scatter points
        matched_col = pg.mkColor('#008000')
        unmatched_col = pg.mkColor('#FFA500')

        # Drop the rows with NaN values in the specified columns (i.e. the international catalogue did not have
        # magtnitudes and/or depths)
        matched_df_drop_mdna = self.matched_df.dropna(subset=['mag_match', 'depth_match'], inplace=False)

        plot1 = self.graph_view.addPlot(0, 0, title="Magnitude Difference",
                                             axisItems={'bottom': DateAxisItem(orientation='bottom',
                                                                               utcOffset=0)})
        plot2 = self.graph_view.addPlot(0, 1, title="Depth Difference",
                                             axisItems={'bottom': DateAxisItem(orientation='bottom',
                                                                               utcOffset=0)})
        plot3 = self.graph_view.addPlot(0, 2, title="Magnitude of Unmatched Events",
                                             axisItems={'bottom': DateAxisItem(orientation='bottom',
                                                                               utcOffset=0)})

        self.lastClicked = []

        self.mag_diff_scatter_plot = pg.ScatterPlotItem(pxMode=True)
        self.mag_diff_scatter_plot.sigClicked.connect(self.scatter_point_clicked)
        self.mag_diff_scatter_plot.addPoints(matched_df_drop_mdna['qtime'], matched_df_drop_mdna['mag_diff'], size=9,
                                             brush=matched_col)
        plot1.addItem(self.mag_diff_scatter_plot)

        self.depth_diff_scatter_plot = pg.ScatterPlotItem(pxMode=True)
        self.depth_diff_scatter_plot.sigClicked.connect(self.scatter_point_clicked)
        self.depth_diff_scatter_plot.addPoints(matched_df_drop_mdna['qtime'], matched_df_drop_mdna['depth_diff'],
                                               size=9, brush=matched_col)
        plot2.addItem(self.depth_diff_scatter_plot)

        self.unmatched_mag_scatter_plot = pg.ScatterPlotItem(pxMode=True)
        self.unmatched_mag_scatter_plot.sigClicked.connect(self.scatter_point_clicked)
        self.unmatched_mag_scatter_plot.addPoints(self.oth_not_matched_df['qtime'], self.oth_not_matched_df['mag'],
                                                  size=9, brush=unmatched_col)
        plot3.addItem(self.unmatched_mag_scatter_plot)

    def scatter_point_deselect(self):
        try:
            for p in self.lastClicked:
                p.resetPen()
        except AttributeError:
            pass

    def scatter_point_clicked(self, plot, points):
        for p in self.lastClicked:
            p.resetPen()
        for p in points:
            p.setPen('r', width=2)
        self.lastClicked = points

        # Highlight the event on the map view and table view when scatter point is clicked
        if plot == self.mag_diff_scatter_plot:
            select_row = self.matched_df.loc[self.matched_df['qtime'] == points[0].pos()[0]]
            self.table_view_highlight(self.tbld.matched_event_table_view, select_row.index.tolist()[0])
        elif plot == self.depth_diff_scatter_plot:
            select_row = self.matched_df.loc[self.matched_df['qtime'] == points[0].pos()[0]]
            self.table_view_highlight(self.tbld.matched_event_table_view, select_row.index.tolist()[0])
        elif plot == self.unmatched_mag_scatter_plot:
            select_row = self.oth_not_matched_df.loc[self.oth_not_matched_df['qtime'] == points[0].pos()[0]]
            self.table_view_highlight(self.tbld.oth_event_table_view, select_row.index.tolist()[0])

    def panMap(self, lng, lat):
        # Pan to the double cliked event (from the table child window)
        self.web_view.page().mainFrame().evaluateJavaScript('map.panTo(L.latLng({}, {}));'.format(lat, lng))
        self.web_view.page().mainFrame().evaluateJavaScript('map.setZoom({});'.format(10))

    @QtCore.pyqtSlot(float, float, str, str, int)
    def onMap_marker_selected(self, lat, lng, event_id, df_id, row_index):
        self.table_view_highlight(self.tbl_view_dict[str(df_id)], row_index)


if __name__ == '__main__':
    # proxy = raw_input("Proxy:")
    # port = raw_input("Proxy Port:")
    # networkProxy = QtNetwork.QNetworkProxy(QtNetwork.QNetworkProxy.HttpProxy, proxy, int(port))
    # QtNetwork.QNetworkProxy.setApplicationProxy(networkProxy)

    app = QtGui.QApplication([])
    app.setStyleSheet(qdarkstyle.load_stylesheet(pyside=False))
    w = MainWindow()
    w.raise_()
    app.exec_()
