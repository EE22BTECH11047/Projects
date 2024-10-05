import React, { Component } from 'react';
import NumberFormat from 'react-number-format';

import Storage from '../../components/Storage';
import Chart from '../../components/Chart';

import CONFIG from '../../config';

class Retirement extends Component {
  constructor(props) {
    super(props);

    let expenses_sum = 0;
    let monthly_revenue = 0;
    let monthly_passive_revenue = 0;
    let estimated_worth = 0;
    let generating_assets = 0;

    props.items.map(function (item) {
      let principal = 0;

      if (item.type === CONFIG.ITEM_TYPE.INCOME) {
        principal = ((item.amount * 12) / CONFIG.INCOME_ESTIMATE_APR) * 100;

        monthly_revenue += item.amount;
      } else if (item.type === CONFIG.ITEM_TYPE.ASSET) {
        principal = item.amount;

        estimated_worth += principal;

        if (item.apr > 0) {
          generating_assets += principal;
        }

        monthly_revenue += (principal * item.apr) / 100 / 12;
        monthly_passive_revenue += (principal * item.apr) / 100 / 12;
      } else if (item.type === CONFIG.ITEM_TYPE.EXPENSE) {
        expenses_sum += item.amount;
      }

      return item;
    });

    this.state = {
      all_items: Storage.getItems(null, props.items),
      assets: Storage.getItems(CONFIG.ITEM_TYPE.ASSET, props.items),
      liabilities: Storage.getItems(CONFIG.ITEM_TYPE.LIABILITY, props.items),
      income: Storage.getItems(CONFIG.ITEM_TYPE.INCOME, props.items),
      expenses: Storage.getItems(CONFIG.ITEM_TYPE.EXPENSE, props.items),
      monthly_revenue: monthly_revenue,
      monthly_passive_revenue: monthly_passive_revenue,
      expenses_sum: expenses_sum,
      estimated_worth: estimated_worth,
      generating_assets: generating_assets,
      chart_data: [],
      retirement_years: 0,
      apr: this.getTotalApr(),
      refresh: false,
    };

    this.handleInputChange = this.handleInputChange.bind(this);
  }

  handleInputChange(event) {
    const target = event.target;
    const value = target.type === 'checkbox' ? target.checked : target.value;
    const name = target.name;

    this.setState({
      [name]: value,
      refresh: !this.state.refresh,
    });

    this.setRetirementYears();
  }

  quantile(arr, q) {
    const asc = (arr) => arr.sort((a, b) => a - b);

    const sorted = asc(arr);
    const pos = (sorted.length - 1) * q;
    const base = Math.floor(pos);
    const rest = pos - base;

    if (sorted[base + 1] !== undefined) {
      return sorted[base] + rest * (sorted[base + 1] - sorted[base]);
    } else {
      return sorted[base];
    }
  }

  componentDidMount() {
    this.setRetirementYears();
  }

  getTotalApr() {
    let total_apr = 0;
    let generating_assets_count = 0;
    let arr = [];

    this.props.items.forEach(function (asset) {
      if (asset.apr > 0) {
        generating_assets_count++;

        arr.push(asset.apr);
      }

      total_apr += asset.apr;
    });

    //const q75 = this.quantile(arr, .90);
    //console.log('q75', q75);

    total_apr /= generating_assets_count;
    total_apr = total_apr || 0;

    return total_apr.toFixed(1);
  }

  setRetirementYears() {
    let year = 1;
    let yearly_contribution =
      (this.state.monthly_revenue - this.state.expenses_sum) * 12;
    let yearly_expenses = this.state.expenses_sum * 12;
    let yearly_apr = this.state.apr;
    let yearly_passive_revenue = this.state.monthly_passive_revenue * 12;
    let estimated_worth = this.state.estimated_worth;
    let generating_assets = this.state.generating_assets;
    let data1 = [yearly_passive_revenue];
    let data2 = [estimated_worth];
    let data3 = [generating_assets];

    while (yearly_passive_revenue < yearly_expenses) {
      const appreciation = (yearly_contribution * year * yearly_apr) / 100;

      yearly_passive_revenue += appreciation;
      estimated_worth += yearly_contribution;
      generating_assets += yearly_contribution;

      //console.log(year, 'estimated_worth: ', estimated_worth, 'yearly_contribution: ', yearly_contribution);

      data1.push(yearly_passive_revenue);
      data2.push(estimated_worth);
      data3.push(generating_assets);

      year++;

      if (year > 100) {
        break;
      }
    }

    this.setState({
      chart_data: [
        {
          name: 'Net Worth',
          data: data2,
        },
        {
          name: 'Generating Assets',
          data: data3,
        },
        {
          name: 'Passive Income',
          data: data1,
        },
      ],
      retirement_years: year,
    });
  }

  renderChart() {
    return (
      <React.Fragment>
        {/*<p>Monthly contribution: ${parseInt(this.state.monthly_revenue - this.state.expenses_sum)}</p>*/}
        <center>
          <p>
            Estimated annual return (APR):&nbsp;
            <input
              type="number"
              placeholder="APR"
              name="apr"
              value={this.state.apr}
              min="0"
              max="100"
              onChange={this.handleInputChange}
            />
            %
          </p>
        </center>
        <br />
        <center>
          <h2>
            You can retire in{' '}
            <strong>
              {this.state.retirement_years + new Date().getFullYear()}
            </strong>{' '}
            <small>({this.state.retirement_years} years)</small>
          </h2>
        </center>
        <br />
        {this.state.chart_data.length && (
          <Chart data={this.state.chart_data} refresh={this.state.refresh} />
        )}
      </React.Fragment>
    );
  }

  render() {
    return (
      <React.Fragment>
        <div className="pricing-header px-3 py-3 pt-md-5 pb-md-4 mx-auto text-center">
          <h1 className="display-4">Retirement calculator</h1>
          <p className="lead">
            How far into the future can you retire based on your current assets,
            passive income, and expenses.
          </p>
        </div>

        <div className="container">
          <div className="card-deck mb-3 text-center">
            <div className="card mb-4 shadow-sm">
              <div className="card-header">
                <h4 className="my-0 font-weight-normal">Assets</h4>
              </div>
              <div className="card-body">
                <h1 className="card-title">
                  <NumberFormat
                    value={this.state.estimated_worth}
                    displayType={'text'}
                    thousandSeparator={true}
                    prefix={'$'}
                    decimalScale={2}
                  />
                </h1>
              </div>
            </div>
            <div className="card mb-4 shadow-sm">
              <div className="card-header">
                <h4 className="my-0 font-weight-normal">Passive Income</h4>
              </div>
              <div className="card-body">
                <h1 className="card-title">
                  <NumberFormat
                    value={this.state.monthly_passive_revenue}
                    displayType={'text'}
                    thousandSeparator={true}
                    prefix={'$'}
                    decimalScale={0}
                  />
                  <small className="text-muted">/ mo</small>
                </h1>
              </div>
            </div>
            <div className="card mb-4 shadow-sm">
              <div className="card-header">
                <h4 className="my-0 font-weight-normal">Expenses</h4>
              </div>
              <div className="card-body">
                <h1 className="card-title">
                  <NumberFormat
                    value={this.state.expenses_sum}
                    displayType={'text'}
                    thousandSeparator={true}
                    prefix={'$'}
                    decimalScale={0}
                  />
                  <small className="text-muted">/ mo</small>
                </h1>
              </div>
            </div>
          </div>

          {this.state.monthly_revenue > this.state.expenses_sum &&
            this.renderChart()}
          {this.state.monthly_revenue <= this.state.expenses_sum && (
            <center>
              <h4>
                Income is lower than expenses. Please add income or remove
                expenses.
              </h4>
            </center>
          )}
        </div>
      </React.Fragment>
    );
  }
}

export default Retirement;
