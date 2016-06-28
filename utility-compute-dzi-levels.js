function getMaximumLevel( width, height ) 
{
  return Math.ceil( Math.log( Math.max( width, height ))/Math.LN2 )
}

function computeLevels( width , height, tileSize)
{
  var maxLevel = getMaximumLevel( width, height )
  var columns;
  var rows;

  for( var level = maxLevel; level >= 0; level-- )
  {
    // compute number of rows & columns
    columns = Math.ceil( width / tileSize )
    rows = Math.ceil( height / tileSize )

    document.write( "<br />level " + level + " is " + width + " x " + height
           + " (" + columns + " columns, " + rows + " rows) = "  + columns * rows + " tiles")

    // compute dimensions of next level
    width  = Math.ceil( width / 2 )
    height = Math.ceil( height / 2 )
  }
}
document.write(computeLevels (50972, 20000, 256)); 