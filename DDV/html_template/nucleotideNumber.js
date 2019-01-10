var PRECISION = 10;      // number of decimal places
var viewer = null;
var pointerStatus = "";
var cursor_in_a_title = false;
var ColumnNumber = 0;
var ColumnRemainder = "";
var PositionInColumn = "";
var originalAspectRatio = originalImageHeight / originalImageWidth;
var Nucleotide = "";
var NucleotideY = "-";
var nucNumX = 0;
var nucNumY = 0;

/**
 fasta_sources, each_layout, ContigSpacingJSON, and contigs are a complete set:
 each is a list with one entry per fasta_source file.  Their indices all match and
 can be used in tandem to fetch different piece of information.
* fasta_sources lists the file names also uses as /chunks/ directories
* each_layout has one LayoutLevels array per file that describes the coordinate frame and origin
* ContigSpacingJSON is the individual contig names placed inside that coordinate frame
* Contigs has a separate object for each file in case some files have duplicate contig names
   This is dynamically loaded from getSequence() and stores actual sequences by name
 */
var contigs = fasta_sources.map(function() { return {} });
var visible_seq_obj;
var theSequenceSplit = []; // used globally by density service
var theSequence = "";
var fragmentid = "";
var sequence_data_loaded = 0;
var sequence_data_viewer_initialized = false;
var file_transfer_in_progress = false;

function init_all(){
    /** Iterates through each chromosome container and initializes and OpenSeaDragon
     * view using the source directory specified in 'data-chr-source' */
    $(".chromosome-container").each(function(index, element){
        init($(element).attr("id"), $(element).attr("data-chr-source"));
    });
}

function init(container_id, source_folder) {
    var source = source_folder == "" ? "" : source_folder + "/"; //ensure directories end with a slash
    source += "GeneratedImages/dzc_output.xml";
    viewer = OpenSeadragon({
        id: container_id,
        prefixUrl: "img/",
        showNavigator: true,
        tileSources: [source],
        maxZoomPixelRatio: 20
    });
    viewer.scalebar({
        type: OpenSeadragon.ScalebarType.MAP,
        pixelsPerMeter: 1,
        minWidth: "70px",
        location: OpenSeadragon.ScalebarLocation.BOTTOM_LEFT,
        xOffset: 5,
        yOffset: 10,
        stayInsideImage: false,
        color: "rgb(30, 30, 30)",
        fontColor: "rgb(10, 10, 10)",
        backgroundColor: "rgba(255, 255, 255, 0.5)",
        fontSize: "normal",
        barThickness: 1,
        sizeAndTextRenderer: OpenSeadragon.ScalebarSizeAndTextRenderer.BASEPAIR_LENGTH
    });

    OpenSeadragon.addEvent(viewer.element, "mousemove", function(event){showNucleotideNumber(event, viewer);});

    //copy content of pointed at sequence fragment to result log
    $("body").keyup(function (event) {
        if (theSequence) {
            if (event.keyCode == 88) {
                $("#outfile").prepend("<div class='sequenceFragment'><div style='background-color:#f0f0f0;'>" + fragmentid + "</div>" + theSequence + "</div>");
            }
        }
    });

    $("#SequenceFragmentInstruction").hide();
    $("#getSequenceButton").hide();
}

function numberWithCommas(x) {
    return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
}

function nucleotide_coordinates_to_sequence_index(index_from_xy, source_index){
    cursor_in_a_title = false;
    var contig_name = "";
    var contig_index = "";
    var index_inside_contig = 0;
    var file_coordinates = "";
    for (var i = 0; i < ContigSpacingJSON[source_index].length; i++) {
        var contig = ContigSpacingJSON[source_index][i];
        if (contig.xy_seq_end > index_from_xy) { // we're in range of the right contig
            contig_index = i;
            if (contig.xy_title_start > index_from_xy) { //we overshot and haven't reached title
                break;
            }
            if (contig.xy_seq_start <= index_from_xy) {// cursor is in nucleotide body
                index_inside_contig = index_from_xy - contig.xy_seq_start + 1;
                file_coordinates = contig.nuc_seq_start + index_inside_contig;
                contig_name = contig.name;
                break;
            } else {// cursor is in label
                cursor_in_a_title = true;
                contig_name = contig.name;
                break;
            }
        }
    }
    return {contig_name: contig_name,
        index_inside_contig: index_inside_contig,
        file_coordinates: file_coordinates,
        contig_index: contig_index,
        fasta_index: source_index};//so we can getSequence() the right file
}

/** Mouse cursor logic for the Ideogram peano layout needs to account for fractal
 * direction switching at each layout level.  Otherwise, similar to
 * tiled_layout_mouse_position()
 * Javascript uses gray code logic to track reversing coordinate frames on odd numbers.*/
function peano_mouse_position(nucNumX, nucNumY, layout_levels, source_index) {
    var index_from_xy = 0;
    var xy_remaining = [nucNumX, nucNumY];
    var axis_flipped = [false, false]
    for (var i = layout_levels.length - 1; i >= 0; i--) {
        var level = layout_levels[i];
        var part = i % 2;
        var number_of_full_increments = Math.floor(xy_remaining[part] / level.thickness);
        var partial_is_reversed = axis_flipped[part];

        var relative_progress = partial_is_reversed ? level.modulo - number_of_full_increments - 1 : number_of_full_increments
        // add total nucleotide size for every full increment of this level e.g. Tile Y height
        index_from_xy += level.chunk_size * relative_progress;

        var step_size_correction = number_of_full_increments * level.thickness;
        //subtract the credited coordinates to shift to relative coordinates in that level
        xy_remaining[part] -= step_size_correction;
        //if the number of full_increments is odd, then we landed in a reversed partial segment
        //at this same level
        var this_level_flipped = number_of_full_increments % 2 == 1;
        //Flip the next axis calculation using XOR to stack recursive layers of flipping
        var next_axis = axis_flipped[(part + 1) %2];
        axis_flipped[(part + 1) %2] = this_level_flipped? !next_axis : next_axis; // XOR

        if (xy_remaining[part] >= level.thickness - level.padding && xy_remaining[part] < level.thickness) {
            return "";//check for invalid coordinate (margins)
        }
    }
    var position_info = nucleotide_coordinates_to_sequence_index(index_from_xy, source_index);
    return position_info;
}


function tiled_layout_mouse_position(nucNumX, nucNumY, layout_levels, source_index) {
    //global variable each_layout set by index.html and python generate_html()
    var index_from_xy = 0;
    var xy_remaining = [nucNumX, nucNumY];
    for (var i = layout_levels.length - 1; i >= 0; i--) {
        var level = layout_levels[i];
        var part = i % 2;
        var number_of_full_increments = Math.floor(xy_remaining[part] / level.thickness);
        // add total nucleotide size for every full increment of this level e.g. Tile Y height
        index_from_xy += level.chunk_size * number_of_full_increments;
        //subtract the credited coordinates to shift to relative coordinates in that level
        xy_remaining[part] -= number_of_full_increments * level.thickness;

        if (xy_remaining[part] >= level.thickness - level.padding && xy_remaining[part] < level.thickness) {
            return "";//check for invalid coordinate (margins)
        }
    }
    var position_info = nucleotide_coordinates_to_sequence_index(index_from_xy, source_index);
    return position_info;
}

function showNucleotideNumber(event, viewer) {
    /** getMousePosition() returns position relative to page,
     * while we want the position relative to the viewer
     * element. so subtract the difference.*/
    var pixel = OpenSeadragon.getMousePosition(event).minus(OpenSeadragon.getElementPosition(viewer.element));
    if (!viewer.isOpen()) {
        return;
    }
    var point = viewer.viewport.pointFromPixel(pixel);
    var position_info = {};
    var information_to_show = false;
    cursor_in_a_title = false;

    if ((point.x < 0) || (point.x > 1)) {
        nucNumX = "-";
        Nucleotide = "";
        pointerStatus = "Outside of Image (X)";
    }
    else {
        nucNumX = Math.round(point.x * originalImageWidth - 0.5) - image_origin[0];
    }

    if ((point.y < 0) || (point.y > originalAspectRatio)) {
        nucNumY = "-";
        Nucleotide = "";
        pointerStatus = "Outside of Image (Y)";
    }
    else {
        nucNumY = Math.round(point.y * originalImageWidth - 0.5) - image_origin[1];
    }

    if ((nucNumX != "-") && (nucNumY != "-")) {
        for (var i = 0; i < each_layout.length; i++) {
            var relX = nucNumX - each_layout[i].origin[0]
            var relY = nucNumY - each_layout[i].origin[1]
            if (relX > -1 && relY > -1){
                if(layout_algorithm == 0){
                    position_info = tiled_layout_mouse_position(relX, relY, each_layout[i].levels, i);
                }else{
                    position_info = peano_mouse_position(relX, relY, each_layout[i].levels, i);
                }
                information_to_show = $.isNumeric(position_info.file_coordinates)
                if(information_to_show){
                    break;
                }
            }
        }
    }


    if(cursor_in_a_title){
        document.getElementById("Nucleotide").innerHTML = position_info.contig_name;
    }else{
        var display_number = information_to_show ? position_info.index_inside_contig : "-";
        document.getElementById("Nucleotide").innerHTML = numberWithCommas(display_number);
        var display_file = information_to_show ? fasta_sources[position_info.fasta_index] : "Sequence under Cursor";
        document.getElementById("FileUnderCursor").innerHTML = display_file;
    }
    //show sequence fragment
    if (sequence_data_viewer_initialized) {
        var lineNumber = "";
        if (information_to_show && position_info.index_inside_contig) {
            var columnWidthInNucleotides = each_layout[position_info.fasta_index].levels[0].modulo
            if(layout_algorithm == 1){
                columnWidthInNucleotides = 100;  // override the tiny display (probably 3)
            }
            Nucleotide = position_info.index_inside_contig;
            lineNumber = Math.floor(Nucleotide / columnWidthInNucleotides);
            var remainder = Nucleotide % columnWidthInNucleotides + columnWidthInNucleotides;
            var start = Math.max(0, (lineNumber - 1) * columnWidthInNucleotides); // not before begin of seq
            var stop = Math.max(start, (lineNumber + 2) * columnWidthInNucleotides); //+2 = +1 start then + width of column
            if(lineNumber == 0){ // first line of the contig
                remainder -= columnWidthInNucleotides;
            }
            if(cursor_in_a_title){
                start = Nucleotide - 1;
                stop = Nucleotide;
            }
            if(contigs[position_info.fasta_index].hasOwnProperty(position_info.contig_name)){
                theSequence = contigs[position_info.fasta_index][position_info.contig_name].substring(start, stop);
                //theSequence = theSequence.replace(/\s+/g, '')
                fragmentid = position_info.contig_name + ": (" +
                  numberWithCommas(start + 1) + " - " + numberWithCommas(stop) + ")";
                //#62 BioJS sequence display dynamically changes to match the number of columns in the current layout
                visible_seq_obj.setNumCols(columnWidthInNucleotides);
                visible_seq_obj.setSequence(theSequence, fragmentid);
                visible_seq_obj.setSelection(remainder, remainder);

                $('#SequenceFragmentInstruction').show();
            }else{
                getSequence(position_info.fasta_index, position_info.contig_index)
            }
        }
        else {
            visible_seq_obj.clearSequence("");
            theSequence = "";
            fragmentid = "";
            $('#SequenceFragmentInstruction').hide();
        }
    }

}

function addLoadEvent(func) {
    var oldOnLoad = window.onload;
    if (typeof window.onload != 'function') {
        window.onload = func;
    }
    else {
        window.onload = function () {
            if (oldOnLoad) {
                oldOnLoad();
            }
            func();
        }
    }
}

function loading_function()
{
    var xhr = new window.XMLHttpRequest();
    //Download progress
    xhr.addEventListener("load", function (evt) {
        $("#status").html("Sequence data loaded.  Display of sequence fragments activated.");
        sequence_data_loaded = 1;
    }, false);
    xhr.addEventListener("progress", function (evt) {
        if (evt.lengthComputable) {
            var percentComplete = (evt.loaded / evt.total) * 100;
            //Do something with download progress
            if (percentComplete < 100) {
                $("#status").html("<img src='img/loading.gif' /> Loading sequence data.  " +
                  "Hover your mouse over the sequence you wish to see: " +
                  parseFloat(percentComplete).toFixed(2) + "% complete");
            }
        }
        else {
            $("#status").html("<img src='img/loading.gif' />Loading sequence data  ... [ " +
              parseFloat(evt.loaded / 1048576).toFixed(2) + " MB loaded ]");
        }
    }, false);
    return xhr;
}

function get_all_sequences() {
    var fasta_path = fasta_sources[0];
    for(let [index, contig] of ContigSpacingJSON.entries()){
        getSequence(0, index); // dispatch one request for each contig
    }
}

function getSequence(fasta_index, contig_index) {
    var fasta_path = "chunks/" + fasta_sources[fasta_index] + "/" + contig_index + ".fa";
    if(!file_transfer_in_progress){
        file_transfer_in_progress = true;
        $.ajax({xhr: loading_function,
            type: "GET",
            url: fasta_path,
            contentType: "text/html",
            success: function (sequence_received) {
                file_transfer_in_progress = false;
                read_contigs(sequence_received, fasta_index);
            },
            error: processInitSequenceError
        });
    }
}

function read_contigs(sequence_received, fasta_index) {
    //read_contigs equiv in javascript
    theSequenceSplit = sequence_received.split(/\r?\n(?=>)/);// begin line, caret  ">");
    for (let contig_s of theSequenceSplit) {
        var lines = contig_s.split(/\r?\n/);
        var title = lines[0].slice(1)
        var seq = lines.slice(1).join('');
        contigs[fasta_index][title] = seq;
    }
    return contigs
}
function init_sequence_view() {
    visible_seq_obj = new Biojs.Sequence({
        sequence: "",
        target: "SeqDisplayTarget",
        format: 'FASTA',
        columns: {size: 100, spacedEach: 0},
        formatSelectorVisible: false,
        fontSize: '18px',
    });
    sequence_data_viewer_initialized = true;
    visible_seq_obj.clearSequence("");
    $('#SequenceFragmentInstruction').hide();
}


function processInitSequenceError() {
    file_transfer_in_progress = false;
};

function outputTable() {
    if (each_layout.length){
       $('#outputContainer').append('<table id="output" style="border: 1px solid #000000;"><tr><th id="FileUnderCursor">Nucleotide Number</th><td id="Nucleotide">-</td></tr></table>    '+
      '<div id="getSequenceButton"><br /><a onclick="get_all_sequences()"> Fetch Sequence </a></div>' +
      '<div id="base"></div><div id="SequenceFragmentFASTA" style="height:200px;">' +
        '<div id="SeqDisplayTarget"></div>' +
        '<div id="SequenceFragmentInstruction" style="display: block;margin-top: -25px;">' +
          'Press "x" key using keyboard to copy this fragment to Result Log</div>' +
        '<div id="status_box">' +
          '<span style="font-weight: bolder;color: darkgrey;font-family: sans-serif;">Status: </span>' +
          '<div id="status"></div>' +
        '</div>' +
      '</div>'
       );
    }
    //provide sequence download for minimal support
    $('#outputContainer').append('<p><a href="sources/">Download Source Data Files</a></p>'
    )
}



function uiLoading (message) {
    //var bufferOutFile=$("#outfile").val();
    $("#status").html(" <img src='img/loading.gif' style='float:left;' /><div style='font-size:11pt;padding-top:10px;padding-bottom:15px;'>"+message+"</div>" );
}

function processError() {
    $("#outfile").prepend("<br />Error.  Connection problem. <div class='resultdivider'></div>");
    $("#status").html("Completed with error." );
}

addLoadEvent(outputTable); // Builds HTML
addLoadEvent(init_all);
addLoadEvent(init_sequence_view);  // could be dependent on a button press
// addLoadEvent(get_all_sequences);

