from PyQt4 import QtCore, QtGui, QtWebKit, QtNetwork
import functools
import qdarkstyle
import pyqtgraph as pg
from obspy.clients.fdsn.client import Client
from obspy.core import UTCDateTime
import sys
import pandas as pd
import numpy as np
from DateAxisItem import DateAxisItem


class PandasModel(QtCore.QAbstractTableModel):
    """
    Class to populate a table view with a pandas dataframe
    """
    def __init__(self, data, parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self._data = data

    def rowCount(self, parent=None):
        return len(self._data.values)

    def columnCount(self, parent=None):
        return self._data.columns.size

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if index.isValid():
            if role == QtCore.Qt.DisplayRole:
                return str(self._data.values[index.row()][index.column()])
        return None

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return self._data.columns[col]
        return None


class TableDialog(QtGui.QDialog):
    def __init__(self, parent=None, isc_df=None, oth_df=None, matched_df=None):
        super(TableDialog, self).__init__(parent)

        self.isc_df = isc_df
        self.oth_df = oth_df
        self.matched_df = matched_df

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

        self.isc_model = PandasModel(self.isc_df)
        self.oth_model = PandasModel(self.oth_df)
        self.matched_model = PandasModel(self.matched_df)

        self.isc_event_table_view.setModel(self.isc_model)
        self.oth_event_table_view.setModel(self.oth_model)
        self.matched_event_table_view.setModel(self.matched_model)

        self.setWindowTitle('Tables')
        self.show()


class MainWindow(QtGui.QWidget):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setupUi()
        self.show()
        self.raise_()

        self.read_events()

    def setupUi(self):

        vbox = QtGui.QVBoxLayout()
        self.setLayout(vbox)

        hbox = QtGui.QHBoxLayout()

        plot_events_button = QtGui.QPushButton('Plot Events')
        plotEvents = functools.partial(self.plot_events)
        plot_events_button.released.connect(plotEvents)
        hbox.addWidget(plot_events_button)
        hbox.addStretch(1)

        vbox.addLayout(hbox)

        view = self.view = QtWebKit.QWebView()

        cache = QtNetwork.QNetworkDiskCache()
        cache.setCacheDirectory("cache")
        view.page().networkAccessManager().setCache(cache)
        view.page().networkAccessManager()

        view.page().mainFrame().addToJavaScriptWindowObject("MainWindow", self)
        view.page().setLinkDelegationPolicy(QtWebKit.QWebPage.DelegateAllLinks)
        view.load(QtCore.QUrl('map.html'))
        view.loadFinished.connect(self.onLoadFinished)
        view.linkClicked.connect(QtGui.QDesktopServices.openUrl)

        vbox.addWidget(view)


        self.graph_view = pg.GraphicsLayoutWidget()

        vbox.addWidget(self.graph_view)

    def onLoadFinished(self):
        with open('map.js', 'r') as f:
            frame = self.view.page().mainFrame()
            frame.evaluateJavaScript(f.read())

    @QtCore.pyqtSlot(float, float, str, str, int)
    def onMap_marker_selected(self, lat, lng, event_id, df_id, row_index):
        self.table_view_highlight(self.tbl_view_dict[str(df_id)], row_index)

    def panMap(self, lng, lat):
        self.view.page().mainFrame().evaluateJavaScript('map.panTo(L.latLng({}, {}));'.format(lat, lng))
        self.view.page().mainFrame().evaluateJavaScript('map.setZoom({});'.format(10))

    def update_graph(self):
        self.graph_view.clear()

        # self.plot1 = self.graph_view.addPlot(0, 0, title="Matched Events", axisItems={'bottom': DateAxisItem(orientation='bottom',
        #                                           utcOffset=0)})
        self.plot2 = self.graph_view.addPlot(0, 0, title="Magnitude Difference", axisItems={'bottom': DateAxisItem(orientation='bottom',
                                                  utcOffset=0)})
        self.plot3 = self.graph_view.addPlot(0, 1, title="Depth Difference", axisItems={'bottom': DateAxisItem(orientation='bottom',
                                                  utcOffset=0)})
        self.plot4 = self.graph_view.addPlot(0, 2, title="Magnitude of Unmatched Events", axisItems={'bottom': DateAxisItem(orientation='bottom',
                                                  utcOffset=0)})

        self.mag_diff_scatter_plot = pg.ScatterPlotItem(pxMode=True)
        self.lastClicked = []
        self.mag_diff_scatter_plot.sigClicked.connect(self.scatter_point_clicked)
        self.mag_diff_scatter_plot.addPoints(self.matched_df['qtime'], self.matched_df['mag_diff'], size=8, brush='g')
        self.plot2.addItem(self.mag_diff_scatter_plot)

        self.depth_diff_scatter_plot = pg.ScatterPlotItem(pxMode=True)
        self.lastClicked = []
        self.depth_diff_scatter_plot.sigClicked.connect(self.scatter_point_clicked)
        self.depth_diff_scatter_plot.addPoints(self.matched_df['qtime'], self.matched_df['depth_diff'], size=8, brush='g')
        self.plot3.addItem(self.depth_diff_scatter_plot)

        self.unmatched_mag_scatter_plot = pg.ScatterPlotItem(pxMode=True)
        self.lastClicked = []
        self.unmatched_mag_scatter_plot.sigClicked.connect(self.scatter_point_clicked)
        self.unmatched_mag_scatter_plot.addPoints(self.isc_not_matched_df['qtime'], self.isc_not_matched_df['mag'], size=8, brush='b')
        self.plot4.addItem(self.unmatched_mag_scatter_plot)

    def scatter_point_clicked(self, plot, points):
        for p in self.lastClicked:
            p.resetPen()
        for p in points:
            p.setPen('r', width=2)
        self.lastClicked = points

        if plot == self.mag_diff_scatter_plot:
            self.select_entry = self.matched_df.loc[self.matched_df['qtime'] == points[0].pos()[0]]
            self.table_view_highlight(self.tbld.matched_event_table_view, self.select_entry.index.tolist()[0])
        elif plot == self.depth_diff_scatter_plot:
            self.select_entry = self.matched_df.loc[self.matched_df['qtime'] == points[0].pos()[0]]
            self.table_view_highlight(self.tbld.matched_event_table_view, self.select_entry.index.tolist()[0])
        elif plot == self.unmatched_mag_scatter_plot:
            self.select_entry = self.isc_not_matched_df.loc[self.isc_not_matched_df['qtime'] == points[0].pos()[0]]
            self.table_view_highlight(self.tbld.isc_event_table_view, self.select_entry.index.tolist()[0])

    def plot_events(self):
        #
        # js_call = "addEvent('{event_id}', '{df_id}', {row_index}, {latitude}, {longitude}, '{a_color}', '{p_color}');" \
        #     .format(event_id="event1", df_id="matched", row_index=1, latitude=-30,
        #             longitude=120, a_color="Red",
        #             p_color="Green")
        # self.view.page().mainFrame().evaluateJavaScript(js_call)
        #
        # js_call = "addEvent('{event_id}', '{df_id}', {row_index}, {latitude}, {longitude}, '{a_color}', '{p_color}');" \
        #     .format(event_id="event2", df_id="isc", row_index=2, latitude=-29,
        #             longitude=122, a_color="Red",
        #             p_color="Blue")
        # self.view.page().mainFrame().evaluateJavaScript(js_call)
        #
        # js_call = "addEvent('{event_id}', '{df_id}', {row_index}, {latitude}, {longitude}, '{a_color}', '{p_color}');" \
        #     .format(event_id="event2", df_id="oth", row_index=3, latitude=-27,
        #             longitude=121, a_color="Red",
        #             p_color="Orange")
        # self.view.page().mainFrame().evaluateJavaScript(js_call)
        #
        # js_call = "createOverlay();"
        # self.view.page().mainFrame().evaluateJavaScript(js_call)



        # Plot the matched events
        for row_index, row in self.matched_df.iterrows():
            js_call = "addEvent('{event_id}', '{df_id}', {row_index}, {latitude}, {longitude}, '{a_color}', '{p_color}');" \
                .format(event_id=row['event_id'], df_id="matched", row_index=int(row_index), latitude=row['lat'], longitude=row['lon'], a_color="Red",
                        p_color="Green")
            self.view.page().mainFrame().evaluateJavaScript(js_call)

        # Plot ISC events (not matched)
        for row_index, row in self.isc_not_matched_df.iterrows():
            js_call = "addEvent('{event_id}', '{df_id}', {row_index}, {latitude}, {longitude}, '{a_color}', '{p_color}');" \
                .format(event_id=row['event_id'], df_id="isc", row_index=int(row_index), latitude=row['lat'], longitude=row['lon'], a_color="Red",
                        p_color="Blue")
            self.view.page().mainFrame().evaluateJavaScript(js_call)

        # Plot oth events (not matched)
        for row_index, row in self.oth_not_matched_df.iterrows():
            js_call = "addEvent('{event_id}', '{df_id}', {row_index}, {latitude}, {longitude}, '{a_color}', '{p_color}');" \
                .format(event_id=row['event_id'], df_id="oth", row_index=int(row_index), latitude=row['lat'], longitude=row['lon'], a_color="Red",
                        p_color="Orange")
            self.view.page().mainFrame().evaluateJavaScript(js_call)


        js_call = "createOverlay();"
        self.view.page().mainFrame().evaluateJavaScript(js_call)


        self.tbld = TableDialog(parent=self, isc_df=self.isc_df, oth_df=self.oth_df, matched_df=self.matched_df)

        # Lookup Dictionary for table views
        self.tbl_view_dict = {"isc": self.tbld.isc_event_table_view,
                              "oth": self.tbld.oth_event_table_view,
                              "matched": self.tbld.matched_event_table_view}

        self.tbld.isc_event_table_view.clicked.connect(self.table_view_clicked)
        self.tbld.oth_event_table_view.clicked.connect(self.table_view_clicked)
        self.tbld.matched_event_table_view.clicked.connect(self.table_view_clicked)

        # Double click and zoom is making markers dissapper when clicing on them
        # self.tbld.isc_event_table_view.doubleClicked.connect(self.table_view_dblclicked)
        # self.tbld.oth_event_table_view.doubleClicked.connect(self.table_view_dblclicked)
        # self.tbld.matched_event_table_view.doubleClicked.connect(self.table_view_dblclicked)

        # plot on graph view
        self.update_graph()

    def table_view_highlight(self, focus_widget, row_index):

        if focus_widget == self.tbld.matched_event_table_view:
            self.tbld.matched_event_table_view.selectRow(row_index)

            # get the row of the the ISC and oth matched quake
            isc_index = int(self.matched_df.loc[row_index, 'isc_ind'])
            oth_index = int(self.matched_df.loc[row_index, 'oth_ind'])
            self.tbld.isc_event_table_view.selectRow(isc_index)
            self.tbld.oth_event_table_view.selectRow(oth_index)
            self.selected_row = self.matched_df.loc[row_index]

            # Highlight the marker on the map
            js_call = "highlightEvent('{event_id}');".format(event_id=self.selected_row['event_id'])
            self.view.page().mainFrame().evaluateJavaScript(js_call)

        elif focus_widget == self.tbld.isc_event_table_view:
            self.tbld.isc_event_table_view.selectRow(row_index)
            self.tbld.oth_event_table_view.clearSelection()
            self.tbld.matched_event_table_view.clearSelection()
            self.selected_row = self.isc_df.loc[row_index]

            # Highlight the marker on the map
            js_call = "highlightEvent('{event_id}');".format(event_id=self.selected_row['event_id'])
            self.view.page().mainFrame().evaluateJavaScript(js_call)

        elif focus_widget == self.tbld.oth_event_table_view:
            self.tbld.oth_event_table_view.selectRow(row_index)
            self.tbld.isc_event_table_view.clearSelection()
            self.tbld.matched_event_table_view.clearSelection()
            self.selected_row = self.oth_df.loc[row_index]

            # Highlight the marker on the map
            js_call = "highlightEvent('{event_id}');".format(event_id=self.selected_row['event_id'])
            self.view.page().mainFrame().evaluateJavaScript(js_call)

    def table_view_dblclicked(self):
        self.table_view_clicked()
        self.panMap(self.selected_row['lon'], self.selected_row['lat'])

    def table_view_clicked(self):
        focus_widget = QtGui.QApplication.focusWidget()
        index = focus_widget.selectionModel().selectedRows()[0]

        self.table_view_highlight(focus_widget, index.row())

    def read_events(self):
        isc_catalogue = Client("IRIS")
        oth_catalogue = Client('http://service.scedc.caltech.edu')

        t_start = UTCDateTime('2016-11-17T00:00:00')
        t_end = UTCDateTime('2016-11-18T00:00:00')

        print('\nRequesting Earthquake Catalogues from Remote Servers.....')

        self.isc_cat = isc_catalogue.get_events(starttime=t_start, endtime=t_end)
        self.oth_cat = oth_catalogue.get_events(starttime=t_start, endtime=t_end)

        # create empty data frame
        self.isc_df = pd.DataFrame(data=None, columns=['oth_ind', 'event_id', 'qtime', 'lat', 'lon', 'depth', 'mag'])
        self.oth_df = pd.DataFrame(data=None, columns=['oth_ind', 'event_id', 'qtime', 'lat', 'lon', 'depth', 'mag'])

        # iterate through the events in oth cat
        for _i, event in enumerate(self.oth_cat):
            print "\r     Parsing event from Local Cat", _i + 1, ' of ', len(self.oth_cat), ' ....',
            sys.stdout.flush()

            # Get quake origin info
            origin_info = event.preferred_origin() or event.origins[0]
            mag_info = event.preferred_magnitude() or event.magnitudes[0]

            self.oth_df.loc[_i] = [int(_i), str(event.resource_id.id), int(origin_info.time.timestamp),
                                   origin_info.latitude, origin_info.longitude,
                                   origin_info.depth, mag_info.mag]

        print('\n')
        # iterate through the events in isc cat
        for _i, event in enumerate(self.isc_cat):
            print "\r     Parsing event from ISC Cat", _i + 1, ' of ', len(self.isc_cat), ' ....',
            sys.stdout.flush()
            # Get quake origin info
            origin_info = event.preferred_origin() or event.origins[0]
            mag_info = event.preferred_magnitude() or event.magnitudes[0]
            self.isc_df.loc[_i] = [int(_i), str(event.resource_id.id), int(origin_info.time.timestamp),
                                   origin_info.latitude, origin_info.longitude,
                                   origin_info.depth, mag_info.mag]


        # =====================Finding matching events =======================

        print('\nFinding Matching Events.....')

        def get_isc_match(row):
            temp = self.isc_df_drop.apply(lambda x: abs(x - row), axis=1)
            # NaNs are treated as small
            smallest_temp = temp.nsmallest(1, columns=['lat', 'lon', 'qtime'])

            isc_index = list(smallest_temp.index)[0]

            if smallest_temp['qtime'].item() <= 15 and \
                    (abs(smallest_temp['lon'].item()) <= 1 or np.isnan(smallest_temp['lon'].item())) and \
                    (abs(smallest_temp['lat'].item()) <= 1 or np.isnan(smallest_temp['lat'].item())):
                ret_s = pd.Series([self.isc_df.loc[isc_index]['oth_ind'].item(), self.isc_df.loc[isc_index]['event_id'],
                                   self.isc_df.loc[isc_index]['qtime'].item(), self.isc_df.loc[isc_index]['lat'].item(),
                                   self.isc_df.loc[isc_index]['lon'].item(),
                                   self.isc_df.loc[isc_index]['depth'].item(), self.isc_df.loc[isc_index]['mag'].item(),
                                   smallest_temp['mag'].item(), smallest_temp['depth'].item()],
                                  index=['isc_ind', 'event_id_match', 'qtime_match', 'lat_match', 'lon_match',
                                         'depth_match', 'mag_match', 'mag_diff', 'depth_diff'])
                return (ret_s)

        # Drop the event_id column (strings) from the data frame to apply vectorised function
        self.oth_df_drop = self.oth_df.drop('event_id', axis=1)
        self.isc_df_drop = self.isc_df.drop('event_id', axis=1)
        self.matched_df = pd.concat((self.oth_df, self.oth_df_drop.apply(get_isc_match, axis=1)), axis=1)

        #drop a row from the matched df if isc_ind in matched df is NaN
        # (I.e. there was no matching earthquake in isc cat)
        self.matched_df.dropna(subset=['isc_ind'], inplace=True)
        self.matched_df.reset_index(drop=True, inplace=True)

        # find isc events not matched and oth events not matched
        self.isc_not_matched_df = self.isc_df[~self.isc_df['oth_ind'].isin(self.matched_df['isc_ind'])]
        self.oth_not_matched_df = self.oth_df[~self.oth_df['oth_ind'].isin(self.matched_df['oth_ind'])]


if __name__ == '__main__':
    # proxy = raw_input("Proxy:")
    # port = raw_input("Proxy Port:")
    # networkProxy = QtNetwork.QNetworkProxy(QtNetwork.QNetworkProxy.HttpProxy, proxy, int(port))
    # QtNetwork.QNetworkProxy.setApplicationProxy(networkProxy)

    app = QtGui.QApplication([])
    app.setStyleSheet(qdarkstyle.load_stylesheet(pyside=False))
    w = MainWindow()
    app.exec_()