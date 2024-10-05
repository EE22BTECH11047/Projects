import React, { Component } from 'react';
import Highcharts from 'highcharts';
import highchartsExporting from 'highcharts/modules/exporting';
import HighchartsReact from 'highcharts-react-official';

highchartsExporting(Highcharts);

class Chart extends Component {
  constructor(props) {
    super(props);

    let data = props.data;

    this.state = {
      options: {
        chart: {
          plotBackgroundColor: null,
          plotBorderWidth: null,
          plotShadow: false,
          backgroundColor: props.backgroundColor || 'transparent',
          type: 'area',
        },
        title: {
          text: props.title,
        },
        tooltip: {
          pointFormat: '{series.name} <b>${point.y:,.0f}</b>',
        },
        xAxis: {
          title: {
            text: 'Year',
          },
          allowDecimals: false,
          labels: {
            formatter: function () {
              return this.value; // clean, unformatted number for year
            },
          },
        },
        yAxis: {
          title: {
            text: 'USD',
          },
        },
        plotOptions: {
          area: {
            pointStart: 2020,
            marker: {
              enabled: false,
              symbol: 'circle',
              radius: 2,
              states: {
                hover: {
                  enabled: true,
                },
              },
            },
          },
        },
        series: data,
      },
    };
  }

  componentWillReceiveProps(props) {
    if (props.refresh !== this.props.refresh) {
      this.setState({
        options: {
          series: props.data,
        },
      });
    }
  }

  render() {
    return (
      <HighchartsReact highcharts={Highcharts} options={this.state.options} />
    );
  }
}

export default Chart;
