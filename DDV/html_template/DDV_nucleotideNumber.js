var PRECISION = 10;      // number of decimal places
var viewer = null;
var pointerStatus = "-";
var ColumnNumber = 0;
var ColumnRemainder = "-";
var PositionInColumn = "-";
var ColumnWidthNoPadding= iLineLength * pixelSize;
var iNucleotidesPerColumn = iLineLength * originalImageHeight / pixelSize;
var ColumnWidth = ColumnPadding + ColumnWidthNoPadding;
var originalAspectRatio = originalImageHeight/originalImageWidth;
var Nucleotide = "-";
var NucleotideY = "-";
var nucNumX = 0;
var nucNumY = 0;

var mySequence;
var theSequenceSplit="";
var theSequence="";
var fragmentid="";
var sequence_data_loaded=0;
var sequence_data_viewer_initialized=0;

function init() {
    viewer = OpenSeadragon({
        id: "container",
        prefixUrl: "img/",
        showNavigator: true,
        tileSources: ["GeneratedImages/dzc_output.xml" ],
        maxZoomPixelRatio: 6
    });
    viewer.scalebar({
        type: OpenSeadragon.ScalebarType.MAP,
        pixelsPerMeter:2,
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

    OpenSeadragon.addEvent(viewer.element, "mousemove", showNucleotideNumber);

    //copy content of pointed at sequence fragment to result log
    $('body').keyup(function (event) {
        if (theSequence){
            if (event.keyCode == 88) {
                $("#outfile").prepend("<div class='sequenceFragment'><div style='background-color:#f0f0f0;'>"+fragmentid+"</div>"+theSequence+"</div>");
            }
        }
    });

    $('#SequenceFragmentInstruction').hide();



}

function showNucleotideNumber(event) {


    // getMousePosition() returns position relative to page,
    // while we want the position relative to the viewer
    // element. so subtract the difference.
    var pixel = OpenSeadragon.getMousePosition(event).minus
    (OpenSeadragon.getElementPosition(viewer.element));

    document.getElementById("mousePixels").innerHTML
      = toString(pixel, true);

    if (!viewer.isOpen()) {
        return;
    }

    var point = viewer.viewport.pointFromPixel(pixel);

    document.getElementById("mousePoints").innerHTML
      = toString(point, true);

    document.getElementById("nucleotideNumberX").innerHTML
      = point.x;
    document.getElementById("nucleotideNumberY").innerHTML
      = point.y;


    if ((point.x < 0) || (point.x > 1)) {
        nucNumX="-";
        Nucleotide = "-";
        pointerStatus = "Outside of Image (X)";

    }
    else {
        nucNumX=(point.x * originalImageWidth).toFixed(0);
    }

    if ((point.y < 0) || (point.y > originalAspectRatio)){
        nucNumY="-";
        Nucleotide = "-";
        pointerStatus = "Outside of Image (Y)";
    }
    else {
        nucNumY=(point.y * originalImageWidth).toFixed(0);
    }

    if ((nucNumX != "-")&&(nucNumY != "-")){
        ColumnNumber = Math.floor(nucNumX/ColumnWidth);
        ColumnRemainder = nucNumX % ColumnWidth;

        PositionInColumn = Math.floor(ColumnRemainder / pixelSize) + 1;
        NucleotideY = iLineLength * Math.floor(nucNumY/pixelSize);

        if ((ColumnRemainder <= ColumnWidth) && (ColumnRemainder >= ColumnWidthNoPadding )){
            ColumnNumber = "-";
            Nucleotide="-";
            PositionInColumn="-";
            pointerStatus = "Outside of Image (Inbetween Columns)";
        }
        else {
            Nucleotide = iNucleotidesPerColumn * ColumnNumber + NucleotideY + PositionInColumn;
            if (Nucleotide > ipTotal) {
                //End of Sequence
                Nucleotide = "-";
            }

        }
    }

    document.getElementById("Nucleotide").innerHTML = Nucleotide;

    //show sequence fragment
    if (sequence_data_viewer_initialized){
        var lineNumber="-";
        if ($.isNumeric(Nucleotide)){
            lineNumber=Math.floor (Nucleotide  /iLineLength);
            remainder=Nucleotide % iLineLength;
            if (lineNumber>0){
                theSequence=theSequenceSplit[lineNumber-1]+theSequenceSplit[lineNumber]+theSequenceSplit[lineNumber+1];
                tempTo=((lineNumber+2)*iLineLength);
                if (ipTotal < tempTo){tempTo=ipTotal;}
                fragmentid= "Sequence fragment at ["+Nucleotide+"], showing: ("+((lineNumber-1)*iLineLength+1)+" - "+tempTo+")";
                mySequence.setSequence(theSequence,fragmentid);
                mySequence.setSelection(remainder+iLineLength, remainder+iLineLength);
            }
            else{
                theSequence=theSequenceSplit[lineNumber]+theSequenceSplit[lineNumber+1];
                fragmentid= "Sequence fragment at ["+Nucleotide+"], showing: "+((lineNumber)*iLineLength+1)+" - "+((lineNumber+2)*iLineLength)+")";
                mySequence.setSequence(theSequence,fragmentid);
                mySequence.setSelection(remainder, remainder);
            }

            $('#SequenceFragmentInstruction').show();

        }
        else{
            mySequence.clearSequence("");
            theSequence="";
            fragmentid="";
            $('#SequenceFragmentInstruction').hide();
        }
    }

}

function toString(point, useParens) {
    var x = point.x;
    var y = point.y;

    if (x % 1 || y % 1) {         // if not an integer,
        x = x.toFixed(PRECISION); // then restrict number of
        y = y.toFixed(PRECISION); // decimal places
    }

    if (useParens) {
        return "(" + x + ", " + y + ")";
    } else {
        return x + " x " + y;
    }
}

function addLoadEvent(func) {

    var oldonload = window.onload;
    if (typeof window.onload != 'function') {
        window.onload = func;
    }
    else {
        window.onload = function () {
            if (oldonload) {
                oldonload();
            }
            func();
        }

    }

}


function getSequence() {


    $.ajax({xhr: function()
    {
        var xhr = new window.XMLHttpRequest();
        //Download progress
        xhr.addEventListener("load", function (evt) {
            $("#status").html("Sequence data loaded.  Display of sequence fragments activated.");
            $("#btnCallGCSkew").click(function (event) {
                GenerateGCSkewChart();
            });
            $("#status").append("<div id='gc-skew-plot-button'>Generate GC Skew activated.");
            sequence_data_loaded = 1;
        }, false);
        xhr.addEventListener("progress", function (evt) {
            if (evt.lengthComputable) {
                var percentComplete = (evt.loaded / evt.total) * 100;
                //Do something with download progress
                if (percentComplete < 100) {
                    $("#status").html("<img src='../../loading.gif' /> Loading sequence data: " + parseFloat(percentComplete).toFixed(2) + "% complete");
                }
            }
            else {
                $("#status").html("<img src='../../loading.gif' />Loading sequence data  ... [ " + parseFloat(evt.loaded / 1048576).toFixed(2) + " MB loaded ]");
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

function initSequence (theSequence) {
    theSequenceSplit=theSequence.split("\n");
    theSequenceSplit.splice(0,1);
    mySequence = new Biojs.Sequence({
        sequence : "",
        target : "SequenceFragmentFASTA",
        format : 'FASTA',
        columns : {size:70,spacedEach:0} ,
        formatSelectorVisible: false,
        fontSize: '11px',
    });
    sequence_data_viewer_initialized=1;
    mySequence.clearSequence("");
    $('#SequenceFragmentInstruction').hide();

}

function processInitSequenceError() {
    //do nothing
};

addLoadEvent(init);
addLoadEvent(getSequence);

function outputTable (){
    outputTable1 = '<table id="output" style="border: 1px solid #000000;"><tr><th>Nucleotide Number</th><td id="Nucleotide">-</td></tr></table><div id="base"></div><div id="SequenceFragmentFASTA" style="height:80px;"><div id="SequenceFragmentInstruction">press "x" key using keyboard to copy this fragment to Result Log</div></div>';
    outputTable21 = '<table class="output" style="border: 1px solid #000000;visibility:hidden;display:none;"><tr><th class="name"> </th><th class="value">Pixels</th><th class="value">Points</th></tr>';
    outputTable22 = '<tr><th>Mouse position</th><td id="mousePixels">-</td><td id="mousePoints">-</td></tr><tr><th>X, Y</th><td id="nucleotideNumberX">-</td><td id="nucleotideNumberY">-</td><td></td></tr>';
    outputTable23 = '<tr><th>(X, Y)</th><td id="NucleotideNumberX">-</td><td id="NucleotideNumberY">-</td></tr><tr><th>Column Number</th><td id="ColumnNumber">-</td><td id="ColumnRemainder">-</td></tr>';
    outputTable24 = '<tr><th>Nucleotide Number</th><td id="Nucleotide">-</td><td>-</td></tr><tr><th>Nucleotides in Local Column</th>   <td id="NucleotideY">-</td><td>-</td></tr>';
    outputTable25 = '<tr><th>Position in Column</th><td id="PositionInColumn">-</td><td></td></tr><tr><th>Nucleotides Per Column</th><td id="iNucleotidesPerColumn">-</td><td></td></tr>';
    outputTable26 = '<tr><th>Aspect Ratio</th><td id="aspectRatio">-</td><td></td></tr><tr><th>Viewport dimensions</th><td id="viewportSizePixels">-</td><td id="viewportSizePoints">-</td></tr></table>';
    document.write(outputTable1);
    document.write(outputTable21);
    document.write(outputTable22);
    document.write(outputTable23);
    document.write(outputTable24);
    document.write(outputTable25);
    document.write(outputTable26);
}

function GenerateGCSkewChart() {

    $("#status").html("<img src='../../loading.gif' />Generating GC Skew Plot...");

    $.getScript("../../d3.v3.js", function(){


        sbegin=$("#sbegin").val();
        send=$("#send").val();
        length = send - sbegin;

//10 000 000 > 10 000
        //10 000 000 > 1 000
        //1 000 000 > 1 00
        //100 000 > 50

//set the default gc_skew_window to 10000
        var gc_skew_window =  10000;
        if (length < 100000 ) {gc_skew_window = 50;}
        else if (length < 1000000) {gc_skew_window = 100;}
        else if (length < 10000000) {gc_skew_window = 1000;}

        var step_G=0;
        var step_C=0;
        var step_GC_skew=0;

        $("#outfile").prepend("<div>GC Skew chart [bp "+sbegin+" to "+send+" ]. GC skew window = "+gc_skew_window+"</div><svg id='gcSkewChart' width='800' height='330'></svg>");


        var lineData = jQuery.map( theSequenceSplit, function( item, index ) {

            step_G += (item.match(/G/g) || []).length;
            step_C += (item.match(/C/g) || []).length;

            if (((index*iLineLength) > sbegin) && ((index*iLineLength) < send) && ((index*iLineLength) % gc_skew_window == 0)){

                if ((step_G + step_C)==0){step_GC_skew=0;}
                else {step_GC_skew = (step_G - step_C)/(step_G + step_C);}
                step_G=0;
                step_C=0;
                return ({'x':(index*iLineLength),'y':step_GC_skew});
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
          .style("font-size","12px")
          .text("Position in sequence ");

        vis.append("svg:g")
          .attr("class", "y axis")
          .attr("transform", "translate(" + (MARGINS.left) + ",0)")
          .call(yAxis)
          .append("text")
          .attr("transform", "rotate(-90)")
          .attr("y", 10)
          .style("text-anchor", "end")
          .style("font-size","12px")
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