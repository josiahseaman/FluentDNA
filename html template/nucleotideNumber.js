var PRECISION = 10;      // number of decimal places
var viewer = null;
var pointerStatus = "-";
var ColumnNumber = 0;
var ColumnRemainder = "-";
var PositionInColumn = "-";
var iNucleotidesPerColumn = columnWidthInNucleotides * originalImageHeight;
var ColumnWidth = ColumnPadding + columnWidthInNucleotides;
var originalAspectRatio = originalImageHeight / originalImageWidth;
var Nucleotide = "-";
var NucleotideY = "-";
var nucNumX = 0;
var nucNumY = 0;

var mySequence;
var wholeSequence = "";
var theSequenceSplit = []; // used globally by density service
var theSequence = "";
var fragmentid = "";
var sequence_data_loaded = 0;
var sequence_data_viewer_initialized = 0;

function init_all(){
    /** Iterates through each chromosome container and initializes and OpenSeaDragon
     * view using the source directory specified in 'data-chr-source' */
    $('.chromosome-container').each(function(index, element){
        init($(element).attr('id'), $(element).attr('data-chr-source'))
    });
}

function init(container_id, source_folder) {
    var source = source_folder == "" ? "" : source_folder + "/"; //ensure directories end with a slash
    source += "GeneratedImages/dzc_output.xml";
    var viewer = OpenSeadragon({
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

    OpenSeadragon.addEvent(viewer.element, "mousemove", function(event){showNucleotideNumber(event, viewer)});

    //copy content of pointed at sequence fragment to result log
    $('body').keyup(function (event) {
        if (theSequence) {
            if (event.keyCode == 88) {
                $("#outfile").prepend("<div class='sequenceFragment'><div style='background-color:#f0f0f0;'>" + fragmentid + "</div>" + theSequence + "</div>");
            }
        }
    });

    $('#SequenceFragmentInstruction').hide();
}

function numberWithCommas(x) {
    return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
}

function classic_layout_mouse_position(nucNumX, nucNumY) {
    var Nucleotide = "-";

    ColumnNumber = Math.floor(nucNumX / ColumnWidth);
    ColumnRemainder = nucNumX % ColumnWidth;

    PositionInColumn = ColumnRemainder + 1;
    NucleotideY = columnWidthInNucleotides * nucNumY;

    if ((ColumnRemainder <= ColumnWidth) && (ColumnRemainder >= columnWidthInNucleotides )) {
        ColumnNumber = "-";
        PositionInColumn = "-";
        pointerStatus = "Outside of Image (Inbetween Columns)";
    }
    else {
        Nucleotide = iNucleotidesPerColumn * ColumnNumber + NucleotideY + PositionInColumn;
        if (Nucleotide > ipTotal) {
            //End of Sequence
            Nucleotide = "-";
        }

    }
    return Nucleotide;
}

function tiled_layout_mouse_position(nucNumX, nucNumY) {
    //global variable layout_levels set by Form1.cs
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
            return "-";//check for invalid coordinate (margins)
        }
    }
    return index_from_xy + 1; // nucleotide position is [1] indexed, instead of [0] indexed
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

    if ((point.x < 0) || (point.x > 1)) {
        nucNumX = "-";
        Nucleotide = "-";
        pointerStatus = "Outside of Image (X)";
    }
    else {
        nucNumX = Math.round(point.x * originalImageWidth - 0.5);
    }

    if ((point.y < 0) || (point.y > originalAspectRatio)) {
        nucNumY = "-";
        Nucleotide = "-";
        pointerStatus = "Outside of Image (Y)";
    }
    else {
        nucNumY = Math.round(point.y * originalImageWidth - 0.5);
    }

    if ((nucNumX != "-") && (nucNumY != "-")) {
        if (layoutSelector == 0) {
            Nucleotide = classic_layout_mouse_position(nucNumX, nucNumY);
        } else {
            Nucleotide = tiled_layout_mouse_position(nucNumX, nucNumY);
        }
    }

    document.getElementById("Nucleotide").innerHTML = numberWithCommas(Nucleotide);

    //show sequence fragment
    if (sequence_data_viewer_initialized) {
        var lineNumber = "-";
        if ($.isNumeric(Nucleotide)) {
            lineNumber = Math.floor(Nucleotide / columnWidthInNucleotides);
            var remainder = Nucleotide % columnWidthInNucleotides + columnWidthInNucleotides;
            var start = Math.max(0, (lineNumber - 1) * columnWidthInNucleotides); // not before begin of seq
            var stop = Math.min(ipTotal, (lineNumber + 2) * columnWidthInNucleotides); //+2 = +1 start then + width of column
            theSequence = wholeSequence.substring(start, stop);
            //user visible indices start at 1, not 0
            fragmentid = "Sequence fragment at [" + numberWithCommas(Nucleotide) +
                "], showing: (" + numberWithCommas(start + 1) + " - " + numberWithCommas(stop) + ")";
            mySequence.setSequence(theSequence, fragmentid);
            mySequence.setSelection(remainder, remainder);

            $('#SequenceFragmentInstruction').show();

        }
        else {
            mySequence.clearSequence("");
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

function getSequence() {
    $('#getSequenceButton').hide();

    $.ajax({
        xhr: function () {
            var xhr = new window.XMLHttpRequest();
            //Download progress
            xhr.addEventListener("progress", function (evt) {
                if (evt.lengthComputable) {
                    var percentComplete = (evt.loaded / evt.total) * 100;
                    //Do something with download progress
                    if (percentComplete < 100) {
                        $("#status").html("<img src='" + output_dir + "loading.gif' /> Loading sequence data: " + parseFloat(percentComplete).toFixed(2) + "% complete");
                    }
                    else {
                        $("#status").html("Sequence data loaded.  Display of sequence fragments activated.");
                        $("#btnCallGCSkew").click(function (event) {
                            GenerateGCSkewChart();
                        });
                        $("#status").append("<div id='gc-skew-plot-button'>Generate GC Skew activated.");
                        sequence_data_loaded = 1;
                    }
                }
                else {
                    $("#status").html("<img src='" + output_dir + "loading.gif' />Loading sequence data  ... [ " + parseFloat(evt.loaded / 1048576).toFixed(2) + " MB loaded ]");
                }
            }, false);
            return xhr;
        },
        type: "GET",
        url: direct_data_file,
        contentType: "text/html",
        success: initSequence,
        error: processInitSequenceError
    });
}

function initSequence(theSequence) {
    theSequenceSplit = theSequence.split("\n");//for removing the newlines
    if (theSequenceSplit[0][0] == ">") { //not all files have a header line
        theSequenceSplit.splice(0, 1); //deletes the header line from the fasta file
    }
    wholeSequence = theSequenceSplit.join('');
    mySequence = new Biojs.Sequence({
        sequence: "",
        target: "SequenceFragmentFASTA",
        format: 'FASTA',
        columns: {size: columnWidthInNucleotides, spacedEach: 0},
        formatSelectorVisible: false,
        fontSize: '11px',
    });
    sequence_data_viewer_initialized = 1;
    mySequence.clearSequence("");
    $('#SequenceFragmentInstruction').hide();

}

function processInitSequenceError() {
    //do nothing
};

addLoadEvent(init_all);
if(includeDensity){
    addLoadEvent(getSequence);
}

function outputTable() {
    document.write('<table id="output" style="border: 1px solid #000000;"><tr><th>Nucleotide Number</th><td id="Nucleotide">-</td></tr></table>');
    document.write("<div id='getSequenceButton'><br /><a onclick='getSequence()'> Fetch Sequence </a></div>");
    document.write('<div id="base"></div><div id="SequenceFragmentFASTA" style="height:80px;"><div id="SequenceFragmentInstruction">press "x" key using keyboard to copy this fragment to Result Log</div></div>');
    document.write('<table class="output" style="border: 1px solid #000000;visibility:hidden;display:none;"><tr><th class="name"> </th><th class="value">Pixels</th><th class="value">Points</th></tr>');
    document.write('<tr><th>Mouse position</th><td id="mousePixels">-</td><td id="mousePoints">-</td></tr><tr><th>X, Y</th><td id="nucleotideNumberX">-</td><td id="nucleotideNumberY">-</td><td></td></tr>');
    document.write('<tr><th>(X, Y)</th><td id="NucleotideNumberX">-</td><td id="NucleotideNumberY">-</td></tr><tr><th>Column Number</th><td id="ColumnNumber">-</td><td id="ColumnRemainder">-</td></tr>');
    document.write('<tr><th>Nucleotide Number</th><td id="Nucleotide">-</td><td>-</td></tr><tr><th>Nucleotides in Local Column</th>   <td id="NucleotideY">-</td><td>-</td></tr>');
    document.write('<tr><th>Position in Column</th><td id="PositionInColumn">-</td><td></td></tr><tr><th>Nucleotides Per Column</th><td id="iNucleotidesPerColumn">-</td><td></td></tr>');
    document.write('<tr><th>Aspect Ratio</th><td id="aspectRatio">-</td><td></td></tr><tr><th>Viewport dimensions</th><td id="viewportSizePixels">-</td><td id="viewportSizePoints">-</td></tr></table>');
}

function GenerateGCSkewChart() {

    $("#status").html("<img src='" + output_dir + "loading.gif' />Generating GC Skew Plot...");

    $.getScript(output_dir + "d3.v3.js", function () {


        var sbegin = $("#sbegin").val();
        var send = $("#send").val();
        var sample_length = send - sbegin;

        //10 000 000 > 10 000
        //10 000 000 > 1 000
        //1 000 000 > 1 00
        //100 000 > 50

        //set the default gc_skew_window to 10000
        var gc_skew_window = 10000;
        if (sample_length < 100000) {
            gc_skew_window = 50;
        }
        else if (sample_length < 1000000) {
            gc_skew_window = 100;
        }
        else if (sample_length < 10000000) {
            gc_skew_window = 1000;
        }

        var step_G = 0;
        var step_C = 0;
        var step_GC_skew = 0;

        $("#outfile").prepend("<div>GC Skew chart [bp " + sbegin + " to " + send + " ]. GC skew window = " + gc_skew_window + "</div><svg id='gcSkewChart' width='800' height='330'></svg>");


        var lineData = jQuery.map(theSequenceSplit, function (item, index) {

            step_G += (item.match(/G/g) || []).length;
            step_C += (item.match(/C/g) || []).length;

            if (((index * columnWidthInNucleotides) > sbegin) && ((index * columnWidthInNucleotides) < send) && ((index * columnWidthInNucleotides) % gc_skew_window == 0)) {

                if ((step_G + step_C) == 0) {
                    step_GC_skew = 0;
                }
                else {
                    step_GC_skew = (step_G - step_C) / (step_G + step_C);
                }
                step_G = 0;
                step_C = 0;
                return ({'x': (index * columnWidthInNucleotides), 'y': step_GC_skew});
            }
            else {
                return null;
            }
        });

        var vis = d3.select("#gcSkewChart"),
            WIDTH = 800,
            HEIGHT = 300,
            MARGINS = {
                top: 20,
                right: 20,
                bottom: 20,
                left: 50
            },
            xRange = d3.scale.linear().range([MARGINS.left, WIDTH - MARGINS.right]).domain([d3.min(lineData, function (d) {
                return d.x;
            }),
                d3.max(lineData, function (d) {
                    return d.x;
                })
            ]),

            yRange = d3.scale.linear().range([HEIGHT - MARGINS.top, MARGINS.bottom]).domain([d3.min(lineData, function (d) {
                return d.y;
            }),
                d3.max(lineData, function (d) {
                    return d.y;
                })
            ]),

            xAxis = d3.svg.axis()
                .scale(xRange)
                .tickSize(5)
                .tickSubdivide(true),

            yAxis = d3.svg.axis()
                .scale(yRange)
                .tickSize(5)
                .orient("left")
                .tickSubdivide(true);


        vis.append("svg:g")
            .attr("class", "x axis")
            .attr("transform", "translate(0," + (HEIGHT - MARGINS.bottom) + ")")
            .call(xAxis)
            .append("text")
            .attr("x", 116)
            .attr("y", 40)
            .style("text-anchor", "start")
            .style("font-size", "12px")
            .text("Position in sequence ");

        vis.append("svg:g")
            .attr("class", "y axis")
            .attr("transform", "translate(" + (MARGINS.left) + ",0)")
            .call(yAxis)
            .append("text")
            .attr("transform", "rotate(-90)")
            .attr("y", 10)
            .style("text-anchor", "end")
            .style("font-size", "12px")
            .text("GC-Skew ");

        var lineFunc = d3.svg.line()
            .x(function (d) {
                return xRange(d.x);
            })
            .y(function (d) {
                return yRange(d.y);
            })
            .interpolate('linear');

        vis.append("svg:path")
            .attr("d", lineFunc(lineData))
            .attr("stroke", "blue")
            .attr("stroke-width", 2)
            .attr("fill", "none");

        $("#status").html("GC Skew Plot added to results.");
    });
}