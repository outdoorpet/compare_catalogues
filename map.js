var map = L.map('map').setView([0, 0], 0);
var layer = new L.StamenTileLayer("toner");
map.addLayer(layer);

var matched_group = L.layerGroup();
var isc_group = L.layerGroup();
var oth_group = L.layerGroup();

var events = {};

function addEvent(event_id, df_id, row_index, latitude, longitude, a_color, p_color) {
    var marker = new L.CircleMarker(
        L.latLng(latitude, longitude), {
            radius: 10,
            color: "Black"
    }).on("click", circleClick);

    marker.status = "--";

    marker.myCustomEventID = event_id;
    marker.myCustomDfID = df_id;
    marker.myCustomRowID = row_index;

    map.addLayer(marker);

    events[event_id] = {
        "marker": marker,
        "latitude": latitude,
        "longitude": longitude,
        "active_color": a_color,
        "passive_color": p_color};

    setMarkerInactive(events[event_id]);

    if (df_id == "matched") {
        matched_group.addLayer(marker);
    } else if (df_id == "isc") {
        isc_group.addLayer(marker);
    } else if (df_id == "oth") {
        oth_group.addLayer(marker);
    }

}

function createOverlay() {
    var overlayMaps = {
        "ISC Events": isc_group,
        "Local Events": oth_group,
        "Matched Events": matched_group};
    L.control.layers(null, overlayMaps).addTo(map);
}

function setMarkerActive(value) {
    if (value.marker.status != "active") {
        value.marker.setStyle({color: value.active_color, opacity: 0.7, fillOpacity: 0.4});
        value.marker.bringToFront()
        value.marker.status = "active";
    }
}


function setMarkerInactive(value) {
    if (value.marker.status != "passive") {
        value.marker.setStyle({color: value.passive_color});
        value.marker.status = "passive";
    }
}


function setAllInactive() {
    _.forEach(events, function(value, key) {
        setMarkerInactive(value);
    });
}


function setAllActive() {
    _.forEach(events, function(value, key) {
        setMarkerActive(value);
    });
}


function highlightEvent(event_id) {
    setAllInactive();
    var value = events[event_id];
    setMarkerActive(value)
}

function resetMarkerSize() {
    _.forEach(events, function(value, key) {
        value.marker.setRadius(10);
    });
}

function arrange_marker_by_depth(sizes) {
    var i = 0;
    _.forEach(events, function(value, key) {
        if (value.marker.myCustomDfID == "matched") {
            value.marker.setRadius(sizes[i]);
            i++;
        }
    });
}


if(typeof MainWindow != 'undefined') {
    function circleClick(e) {
        var clickedCircle = e.target;

//    highlightEvent(clickedCircle.myCustomID)
    clickedCircle.bindPopup(clickedCircle.myCustomEventID).openPopup();
    MainWindow.onMap_marker_selected(clickedCircle.getLatLng().lat, clickedCircle.getLatLng().lng, clickedCircle.myCustomEventID, clickedCircle.myCustomDfID, clickedCircle.myCustomRowID)
    }

//    function overlay_select
}