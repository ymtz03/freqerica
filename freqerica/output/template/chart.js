// Correlation Function
( function() {

  //d3.select("body").style("margin", 0);
  
  var width = 800, height = 400;
  var svg = d3.select("#correlationfunction")
              //.style("background-color", "rgb(240, 255, 255)").style("padding-left", 4)
              .append("svg").attr("viewBox", [0,0, width, height]);
  var margin = {top: 20, right:40, bottom:30, left:40};

  /*
  svg.append("rect")
     .attr("x", 0)
     .attr("y", 0)
     .attr("width", width)
     .attr("height", height)
     .attr("fill","white");

  svg.append("rect")
     .attr("x", margin.left)
     .attr("y", margin.top)
     .attr("width", width-(margin.right + margin.left))
     .attr("height", height-(margin.top+margin.bottom))
     .attr("fill","white");
  */

  var time = data.corr_re_exact.map((d,i)=>i*data.dt);
  console.log(time);

  var x = d3.scaleLinear()
            .domain(d3.extent(time))
            .range([margin.left, width-margin.right])
            .nice();

  var y = d3.scaleLinear()
            .domain([-1,1])
            .range([height-margin.bottom, margin.top])
            .nice();

  var vgrid = svg.append("g")
                 .selectAll("path")
                 .data(x.ticks())
                 .enter()
                 .append("line")
                 .attr("stroke-width", 1)
                 .attr("stroke-dasharray", "4")
                 .attr("stroke", "lightgray")
                 .attr("x1", d=>x(d))
                 .attr("x2", d=>x(d))
                 .attr("y1", height-margin.bottom)
                 .attr("y2", margin.top);

  var line = d3.line().x((d,i)=>x(time[i])).y(d=>y(d));

  svg.append("path")
     .datum(data.corr_re_trott)
     .attr("fill", "none")
     .attr("stroke", "black")
     .attr("stroke-width", 1.5)
     .attr("d", line);
  
  svg.append("path")
     .datum(data.corr_im_trott)
     .attr("fill", "none")
     .attr("stroke", "black")
     .attr("stroke-width", 1.5)
     .attr("stroke-dasharray", '4')
     .attr("d", line);

  svg.append("path")
     .datum(data.corr_re_exact)
     .attr("fill", "none")
     .attr("stroke", "steelblue")
     .attr("stroke-width", 1.5)
     .attr("d", line);

  svg.append("path")
     .datum(data.corr_im_exact)
     .attr("fill", "none")
     .attr("stroke", "steelblue")
     .attr("stroke-width", 1.5)
     .attr("stroke-dasharray", '4')
     .attr("d", line);

  
  var xaxis = svg.append("g")
                 .call(d3.axisBottom(x))
                 .attr("transform", "translate(0," + y(0) + ")" );
                 //.attr("transform", "translate(0," + (height-margin.bottom) + ")" );

  var yaxis = svg.append("g")
                 .call(d3.axisLeft(y))
                 .attr("transform", "translate(" + (margin.left) + ",0)" );

})();

// Spectrum
( function() {

  var width = 800, height = 400;
  var svg = d3.select("#spectrum").append("svg").attr("viewBox", [0,0, width, height]);
  var margin = {top: 20, right:40, bottom:30, left:40};

  
/*
  svg.append("rect")
     .attr("x", 0)
     .attr("y", 0)
     .attr("width", width)
     .attr("height", height)
     .attr("fill","white");

  svg.append("rect")
     .attr("x", margin.left)
     .attr("y", margin.top)
     .attr("width", width-(margin.right + margin.left))
     .attr("height", height-(margin.top+margin.bottom))
     .attr("fill","white");
  */

  var x = d3.scaleLinear()
            .domain(d3.extent(data.energy_range))
            .range([margin.left, width-margin.right])
            .nice();

  var y = d3.scaleLinear()
            .domain(d3.extent(data.spectrum))
            .range([height-margin.bottom, margin.top])
            .nice();

  var y_prony = d3.scaleLinear()
                  .domain([0,1])
                  .range([height-margin.bottom, margin.top])
                  .nice();

  var vgrid = svg.append("g")
                 .selectAll("path")
                 .data(x.ticks())
                 .enter()
                 .append("line")
                 .attr("stroke-width", 1)
                 .attr("stroke-dasharray", "4")
                 .attr("stroke", "lightgray")
                 .attr("x1", d=>x(d))
                 .attr("x2", d=>x(d))
                 .attr("y1", height-margin.bottom)
                 .attr("y2", margin.top);

  var eigvals = svg.append("g")
                   .selectAll("path")
                   .data(data.energy)
                   .enter()
                   .append("line")
                   .attr("stroke-width", 1.5)
                   .attr("stroke", "steelblue")
                   .attr("x1", d=>x(d))
                   .attr("x2", d=>x(d))
                   .attr("y1", height-margin.bottom)
                   .attr("y2", margin.top);

  var prony_exact = svg.append("g").selectAll("xxx").data(data.prony_p_exact).enter().append("circle")
                       .attr("cx", d=>x((-d-2*Math.PI)/data.dt) )
                       .attr("cy", (d,i)=>y_prony(data.prony_A_exact[i]))
                       .attr("r", 3)
                       .attr("stroke","red")
                       .attr("stroke-width",1)
                       .attr("fill","none");

  var prony_trott = svg.append("g").selectAll("xxx").data(data.prony_p_trott).enter().append("circle")
                       .attr("cx", d=>x((-d-2*Math.PI)/data.dt) )
                       .attr("cy", (d,i)=>y_prony(data.prony_A_trott[i]))
                       .attr("r", 3)
                       .attr("stroke","green")
                       .attr("stroke-width",1)
                       .attr("fill","none");

  var freq_vs_amp = [];
  for(var i=0; i<data.energy_range.length; ++i) {
    freq_vs_amp.push([data.energy_range[i], data.spectrum[i]]);
  }

  var line = d3.line().x(d=>x(d[0])).y(d=>y(d[1]));

  svg.append("path")
     .datum(freq_vs_amp)
     .attr("fill", "none")
     .attr("stroke", "black")
     .attr("stroke-width", 1.5)
     .attr("d", line);

  var xaxis = svg.append("g")
                 .call(d3.axisBottom(x))
                 .attr("transform", "translate(0," + (height-margin.bottom) + ")" );

  var yaxis = svg.append("g")
                 .call(d3.axisLeft(y))
                 .attr("transform", "translate(" + (margin.left) + ",0)" );

  var yaxis = svg.append("g")
                 .call(d3.axisRight(y_prony))
                 .attr("transform", "translate(" + (width-margin.right) + ",0)" );
})();
